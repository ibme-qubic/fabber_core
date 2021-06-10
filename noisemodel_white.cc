/*  noisemodel_white.cc - Class implementation for the multiple white noise model

 Adrian Groves and Michael Chappell, FMRIB Image Analysis Group & IBME QuBIc Group

 Copyright (C) 2007-2015 University of Oxford  */

/*  CCOPYRIGHT */

#include "noisemodel_white.h"

#include "easylog.h"
#include "noisemodel.h"
#include "rundata.h"
#include "tools.h"

#include <miscmaths/miscmaths.h>
#include "armawrap/newmat.h"

#include <ostream>
#include <string>

using MISCMATHS::digamma;
using namespace NEWMAT;
using namespace std;

NoiseModel *WhiteNoiseModel::NewInstance()
{
    return new WhiteNoiseModel();
}

WhiteParams::WhiteParams(int N)
    : nPhis(N)
    , phis(N)
{
}

WhiteParams::WhiteParams(const WhiteParams &from)
    : nPhis(from.nPhis)
    , phis(from.phis)
{
}

WhiteParams *WhiteParams::Clone() const
{
    return new WhiteParams(*this);
}
const WhiteParams &WhiteParams::operator=(const NoiseParams &in)
{
    const WhiteParams &from = dynamic_cast<const WhiteParams &>(in);
    assert(nPhis == from.nPhis);
    phis = from.phis;
    return *this;
}

const MVNDist WhiteParams::OutputAsMVN() const
{
    assert((unsigned)nPhis == phis.size());
    MVNDist mvn(phis.size());
    SymmetricMatrix vars(phis.size());
    vars = 0;
    for (unsigned i = 1; i <= phis.size(); i++)
    {
        mvn.means(i) = phis[i - 1].CalcMean();
        vars(i, i) = phis[i - 1].CalcVariance();
    }
    mvn.SetCovariance(vars);
    return mvn;
}

void WhiteParams::InputFromMVN(const MVNDist &mvn)
{
    for (unsigned i = 1; i <= phis.size(); i++)
    {
        phis[i - 1].SetMeanVariance(mvn.means(i), mvn.GetCovariance()(i, i));
        for (int j = i + 1; j <= mvn.means.Nrows(); j++)
            if (mvn.GetCovariance()(i, j) != 0.0)
                throw FabberRunDataError("Phis should have zero covariance!");
    }
}

void WhiteParams::Dump(ostream &os) const
{
    assert((unsigned)nPhis == phis.size());
    for (unsigned i = 0; i < phis.size(); i++)
    {
        os << "WhiteNoiseModel::Phi_" << i + 1 << ": ";
        phis[i].Dump(os);
    }
}

void WhiteNoiseModel::Initialize(FabberRunData &args)
{
    NoiseModel::Initialize(args);

    // White noise can have a pattern. This is represented as, e.g.
    // 123123123... or 123456789ab123456789ab...
    // Whatever the patter is, each distinct digit/letter defines
    // a noise parameter, and the pattern is repeated along the
    // data.
    phiPattern = args.GetStringDefault("noise-pattern", "1");
    assert(phiPattern.length() > 0);

    // This is just a quick way to validate the input. It will not
    // set up the pattern correctly as that requires the data length
    // (number of timeseries samples) to be known and is done
    // at calculation time.
    MakeQis(phiPattern.length());

    // Allow phi to be locked externally
    lockedNoiseStdev = convertTo<double>(args.GetStringDefault("locked-noise-stdev", "-1"));
    assert(lockedNoiseStdev == -1 || lockedNoiseStdev > 0); // TODO (should throw)

    // Set phi prior externally
    phiprior = convertTo<double>(args.GetStringDefault("prior-noise-stddev", "-1"));
    if (phiprior < 0 && phiprior != -1)
        throw InvalidOptionValue("prior-noise-stddev", stringify(phiprior), "Must be > 0");
}

int WhiteNoiseModel::NumParams()
{
    return Qis.size();
}
WhiteParams *WhiteNoiseModel::NewParams() const
{
    return new WhiteParams(Qis.size());
}
void WhiteNoiseModel::HardcodedInitialDists(NoiseParams &priorIn, NoiseParams &posteriorIn) const
{
    WhiteParams &prior = dynamic_cast<WhiteParams &>(priorIn);
    WhiteParams &posterior = dynamic_cast<WhiteParams &>(posteriorIn);

    int nPhis = Qis.size();
    assert(nPhis > 0);
    //    prior.resize(nPhis);
    //    posterior.resize(nPhis);
    assert(nPhis == prior.nPhis);
    assert(nPhis == posterior.nPhis);

    for (int i = 1; i <= nPhis; i++)
    {
        if (phiprior == -1)
        {
            // use non-informative prior on phi
            posterior.phis[i - 1].b = prior.phis[i - 1].b = 1e6;
            posterior.phis[i - 1].c = prior.phis[i - 1].c = 1e-6;

            // Bit of black magic... tiny initial noise precision seems to help
            posterior.phis[i - 1].b = 1e-8;
            posterior.phis[i - 1].c = 50;
        }
        else
        {
            // use input nosie std dev to determine the prior (and inital values) for phi
            // this is for the case where N is small, so it make sense to use a large(ish) c_prior
            // since C ~ N (data points) and we struggle to estimate noise when N is small.
            posterior.phis[i - 1].c = prior.phis[i - 1].c
                = 0.5; // assumes that a given std deviation is equivelent to N=1 measurements
            posterior.phis[i - 1].b = prior.phis[i - 1].b
                = 1 / (phiprior * phiprior
                          * prior.phis[i - 1]
                                .c); // NB phiprior is a std dev for the noise that is read in
        }
    }
}

void WhiteNoiseModel::MakeQis(int dataLen) const
{
    if (!Qis.empty() && Qis[0].Nrows() == dataLen)
        return; // Qis are already up-to-date

    // Read the pattern string into a vector pat
    const int patternLen = phiPattern.length();
    assert(patternLen > 0);
    if (patternLen > dataLen)
        throw InvalidOptionValue("noise-pattern", phiPattern, "Pattern length exceeds data length");

    vector<int> pat;
    int nPhis = 0;

    for (int i = 1; i <= patternLen; i++)
    {
        char c = phiPattern[i - 1];
        int n;
        if (c >= '1' && c <= '9') // Index from 1, so 0 is not allowed
            n = (c - '0');
        else if (c >= 'A' && c <= 'Z')
            n = (c - 'A' + 10);
        else if (c >= 'a' && c <= 'z')
            n = (c - 'a' + 10);
        else
            throw InvalidOptionValue("noise-pattern", stringify(c), "Invalid character");
        pat.push_back(n);
        if (nPhis < n)
            nPhis = n;
    }

    // Extend pat to the the full data length by repeating
    // Inefficient but pretty harmless unless you have huge timeseries
    while ((int)pat.size() < dataLen)
        pat.push_back(pat.at(pat.size() - patternLen));

    LOG << "WhiteNoiseMode::Pattern of phis used is " << pat << endl;

    // Regenerate Qis. This is a vector of
    // diagonal matrices, one for each parameter (Phi).
    // Each matrix is the same size as the length of the data.
    Qis.clear();
    DiagonalMatrix zeroes(dataLen);
    zeroes = 0.0;
    Qis.resize(nPhis, zeroes); // initialize all to zero

    // For each sample in the timeseries, find the
    // appropriate parameter (phi) from the pattern
    // and set the diagonal element in the Qi matrix
    // to 1. So each Qi matrix has elements to define
    // what samples it applies to
    for (int d = 1; d <= dataLen; d++)
    {
        // Only flag a time point as relevant if it is not masked
        if (std::find(m_masked_tpoints.begin(), m_masked_tpoints.end(), d)
            == m_masked_tpoints.end())
        {
            Qis.at(pat.at(d - 1) - 1)(d, d) = 1;
        }
    }
}

void WhiteNoiseModel::UpdateNoise(NoiseParams &noise, const NoiseParams &noisePrior,
    const MVNDist &theta, const LinearFwdModel &linear, const ColumnVector &data) const
{
    WhiteParams &posterior = dynamic_cast<WhiteParams &>(noise);
    const WhiteParams &prior = dynamic_cast<const WhiteParams &>(noisePrior);

    const Matrix &J = linear.Jacobian();
    ColumnVector k = data - linear.Offset() + J * (linear.Centre() - theta.means);

    // Check there are the same number of Qis in this model and in the
    // prior and posterior parameter sets.
    MakeQis(data.Nrows());
    const int nPhis = Qis.size();
    assert(nPhis == posterior.nPhis);
    assert(nPhis == prior.nPhis);

    // Update each phi distribution in turn
    for (int i = 1; i <= nPhis; i++)
    {
        // Each Phi matrix is a diagonal matrix of same size size as the
        // number of time samples in the data
        const DiagonalMatrix &Qi = Qis[i - 1];

        // This is calculating the 2nd and 3rd terms of RHS of Eq (22) in Chappel et al 2009
        double tmp = (k.t() * Qi * k).AsScalar() + (theta.GetCovariance() * J.t() * Qi * J).Trace();

        // This is Eq (22) in Chappel et al 2009
        posterior.phis[i - 1].b = 1 / (tmp * 0.5 + 1 / prior.phis[i - 1].b);

        // Number of data sample points which use this parameter.
        // Should be an integer
        double nTimes = Qi.Trace();
        assert(nTimes == int(nTimes));

        // This is Eq (21) in Chappel et al 2009
        posterior.phis[i - 1].c = (nTimes - 1) * 0.5 + prior.phis[i - 1].c;

        if (lockedNoiseStdev > 0)
        {
            // Ignore this update and force phi to a specified value.
            // b*c = noise precision = lockedNoiseStdev^-2
            posterior.phis[i - 1].b
                = 1 / posterior.phis[i - 1].c / lockedNoiseStdev / lockedNoiseStdev;
        }
    }
}

void WhiteNoiseModel::UpdateTheta(const NoiseParams &noiseIn, MVNDist &theta,
    const MVNDist &thetaPrior, const LinearFwdModel &linear, const ColumnVector &data,
    MVNDist *thetaWithoutPrior, float LMalpha) const
{
    const WhiteParams &noise = dynamic_cast<const WhiteParams &>(noiseIn);

    const ColumnVector &ml = linear.Centre();
    const ColumnVector &gml = linear.Offset();
    const Matrix &J = linear.Jacobian();

    // Make sure Qis are up-to-date
    MakeQis(data.Nrows());
    assert(Qis.size() == (unsigned)noise.nPhis);

    // Marginalize over phi distributions
    // Qis are diagonal matrices with 1 only where that phi applies.
    // Adding up all the Qis will give you the identity matrix.
    DiagonalMatrix X(data.Nrows());
    X = 0;
    for (unsigned i = 1; i <= Qis.size(); i++)
        X += Qis[i - 1] * noise.phis[i - 1].CalcMean();

    // Update Lambda (model precisions)
    //
    // This is Eq (19) in Chappel et al (2009)
    //
    // use << instead of = because this is considered a lossy assignment
    // (since NEWMAT isn't smart enough to know J'*X*J is always symmetric)
    SymmetricMatrix Ltmp;
    Ltmp << J.t() * X * J;
    theta.SetPrecisions(thetaPrior.GetPrecisions() + Ltmp);

    // Error checking
    LogAndSign chk = theta.GetPrecisions().LogDeterminant();
    if (chk.Sign() <= 0)
    {
        LOG << "WhiteNoiseModel:: In UpdateTheta, theta precisions aren't positive-definite: "
            << chk.Sign() << ", " << chk.LogValue() << endl;
        LOG << "Means: " << theta.means.t() << endl;
        LOG << "Precisions: " << endl << theta.GetPrecisions() << endl;
        LOG << "Data: " << data.t() << endl;
    }

    // Update m (model means)
    //
    // This is the first term of RHS of Eq (20) in Chappel et al (2009)
    ColumnVector mTmp = J.t() * X * (data - gml + J * ml);
    if (LMalpha <= 0.0)
    {
        // Normal update (NB the LM update reduces to this when alpha=0 strictly)
        // This is Eq (20) in Chappel et al (2009). Note that covariance of theta
        // is inverse of precisions.
        theta.means
            = theta.GetCovariance() * (mTmp + thetaPrior.GetPrecisions() * thetaPrior.means);
    }
    else
    {
        // We are in LM mode so use the appropriate update
        // See Appendix C in Chappel et al (2009)
        Matrix Delta;
        SymmetricMatrix prec;
        DiagonalMatrix precdiag;
        prec = theta.GetPrecisions();
        precdiag << prec;

        // a different (but equivalent) form for the LM update
        Delta = J.t() * X * (data - gml) + thetaPrior.GetPrecisions() * thetaPrior.means
            - thetaPrior.GetPrecisions() * ml;
        try
        {
            theta.means = ml + (prec + LMalpha * precdiag).i() * Delta;
        }
        catch (Exception)
        {
            WARN_ONCE("WhiteNoiseMode: matrix was singular in LM update");
        }
        // LM update - old method
        // theta.means = (prec + LMalpha*precdiag).i()
        // * ( mTmp + thetaPrior.GetPrecisions() * thetaPrior.means );
    }

    // Optional update of model parameters without covariance rior
    if (thetaWithoutPrior != NULL)
    {
        thetaWithoutPrior->SetSize(theta.GetSize());
        thetaWithoutPrior->SetPrecisions(Ltmp);
        thetaWithoutPrior->means = thetaWithoutPrior->GetCovariance() * mTmp;
    }
}

double WhiteNoiseModel::CalcFreeEnergy(const NoiseParams &noiseIn, const NoiseParams &noisePriorIn,
    const MVNDist &theta, const MVNDist &thetaPrior, const LinearFwdModel &linear,
    const ColumnVector &data) const
{
    const int nPhis = Qis.size();
    const WhiteParams &noise = dynamic_cast<const WhiteParams &>(noiseIn);
    const WhiteParams &noisePrior = dynamic_cast<const WhiteParams &>(noisePriorIn);

    // Calculate some matrices we will need
    const Matrix &J = linear.Jacobian();
    ColumnVector k = data - linear.Offset() + J * (linear.Centre() - theta.means);
    const SymmetricMatrix &Linv = theta.GetCovariance();

    // some values we will need
    int nTimes = data.Nrows()
        - m_masked_tpoints.size(); //*NB assume that each row is an individual time point
    int nTheta = theta.means.Nrows();

    // The following is based on noisemodel_ar::CalcFreeEnergy, modified to remove ar parts - MAC
    // 11-7-2007. Some modifications have been made for consistency with (MAC)varbayes2.m - these
    // are noted

    // calcualte individual parts of the free energy
    double expectedLogThetaDist = // bits arising from the factorised posterior for theta
        +0.5 * theta.GetPrecisions().LogDeterminant().LogValue()
        - 0.5 * nTheta * (log(2 * M_PI) + 1);

    double expectedLogPhiDist = 0; // bits arising fromt he factorised posterior for phi
    vector<double> expectedLogPosteriorParts(10); // bits arising from the likelihood
    for (int i = 0; i < 10; i++)
        expectedLogPosteriorParts[i] = 0;

    for (int i = 0; i < nPhis; i++)
    {
        double si = noise.phis[i].b;
        double ci = noise.phis[i].c;
        double siPrior = noisePrior.phis[i].b;
        double ciPrior = noisePrior.phis[i].c;

        expectedLogPhiDist += -gammaln(ci) - ci * log(si) - ci + (ci - 1) * (digamma(ci) + log(si));

        expectedLogPosteriorParts[0] += (digamma(ci) + log(si))
            * ((Qis[i].Trace()) * 0.5 + ciPrior - 1); // nTimes using phi_{i+1} = Qis[i].Trace()

        expectedLogPosteriorParts[9]
            += -gammaln(ciPrior) - ciPrior * log(siPrior) - si * ci / siPrior;
    }

    expectedLogPosteriorParts[1] = 0; //*NB not required

    expectedLogPosteriorParts[2]
        = -0.5 * (k.t() * k).AsScalar() - 0.5 * (J.t() * J * Linv).Trace(); //*NB remove Qsum

    expectedLogPosteriorParts[3] = +0.5 * thetaPrior.GetPrecisions().LogDeterminant().LogValue()
        - 0.5 * nTimes * log(2 * M_PI) - 0.5 * nTheta * log(2 * M_PI);

    expectedLogPosteriorParts[4] = -0.5
        * ((theta.means - thetaPrior.means).t() * thetaPrior.GetPrecisions()
              * (theta.means - thetaPrior.means))
              .AsScalar();

    expectedLogPosteriorParts[5] = -0.5 * (Linv * thetaPrior.GetPrecisions()).Trace();

    expectedLogPosteriorParts[6] = 0; //*NB not required

    expectedLogPosteriorParts[7] = 0; //*NB not required

    expectedLogPosteriorParts[8] = 0; //*NB not required

    // Assemble the parts into F
    double F = -expectedLogThetaDist - expectedLogPhiDist;

    for (int i = 0; i < 10; i++)
        F += expectedLogPosteriorParts[i];

    // Error checking
    if (!(F - F == 0))
    {
        LOG_ERR("WhiteNoiseModel::expectedLogThetaDist == " << expectedLogThetaDist << endl);
        LOG_ERR("WhiteNoiseModel::expectedLogPhiDist == " << expectedLogPhiDist << endl);
        // LOG_ERR("expectedLogPosteriorParts == " << expectedLogPosteriorParts << endl);
        throw FabberInternalError("WhiteNoiseModel::Non-finite free energy!");
    }

    return F;
}
