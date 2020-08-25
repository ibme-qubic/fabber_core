/**
 * prior.cc
 *
 * Class for a parameter prior
 *
 * Copyright (C) 2007-2017 University of Oxford
 */

/*  CCOPYRIGHT */

#include "priors.h"

#include "dist_mvn.h"
#include "rundata.h"
#include "tools.h"

#include <miscmaths/miscmaths.h>
#include <newmat.h>

#include <math.h>
#include <ostream>
#include <string>
#include <vector>

using namespace std;
using namespace NEWMAT;
using MISCMATHS::digamma;

std::ostream &operator<<(std::ostream &out, const Prior &prior)
{
    prior.DumpInfo(out);
    return out;
}

string Prior::ExpandPriorTypesString(string priors_str, unsigned int num_params)
{
    // Find out how many prior types are in the string, and what the + character
    // should be interpreted as
    unsigned int n_str_params = 0;
    char repeat_type = '-';
    bool plus_found = false;
    for (size_t i = 0; i < priors_str.size(); i++)
    {
        if (priors_str[i] != '+')
        {
            if (!plus_found)
                repeat_type = priors_str[i];
            n_str_params++;
        }
        else if (plus_found)
        {
            throw InvalidOptionValue(
                "param-spatial-priors", priors_str, "Only one + character allowed");
        }
        else
        {
            plus_found = true;
        }
    }

    if (n_str_params > num_params)
    {
        throw InvalidOptionValue("param-spatial-priors", priors_str, "Too many parameters");
    }
    else if (n_str_params < num_params)
    {
        // Expand '+' char, if present, to give correct number of parameters
        // If there is no +, append with '-', meaning 'model default'
        int deficit = num_params - n_str_params;
        size_t plus_pos = priors_str.find("+");
        if (plus_pos != std::string::npos)
        {
            priors_str.insert(plus_pos, deficit - 1, '+');
        }
        else
        {
            priors_str.insert(priors_str.end(), deficit, '-');
        }
    }
    else
    {
        // We already have enough types for all the parameters so erase any
        // pointless + char
        priors_str.erase(std::remove(priors_str.begin(), priors_str.end(), '+'), priors_str.end());
    }

    // Finally, replace all + chars with identified repeat type
    std::replace(priors_str.begin(), priors_str.end(), '+', repeat_type);
    assert(priors_str.size() == num_params);

    return priors_str;
}

DefaultPrior::DefaultPrior(const Parameter &p)
    : m_param_name(p.name)
    , m_idx(p.idx)
    , m_type_code(p.prior_type)
    , m_params(p.prior)
{
}

void DefaultPrior::DumpInfo(std::ostream &out) const
{
    out << "DefaultPrior: Parameter " << m_idx << " '" << m_param_name << "'"
        << " mean: " << m_params.mean() << " precision: " << m_params.prec();
}

double DefaultPrior::ApplyToMVN(MVNDist *prior, const RunContext &ctx)
{
    prior->means(m_idx + 1) = m_params.mean();

    SymmetricMatrix prec = prior->GetPrecisions();
    prec(m_idx + 1, m_idx + 1) = m_params.prec();
    prior->SetPrecisions(prec);

    return 0;
}

ImagePrior::ImagePrior(const Parameter &p, FabberRunData &rundata)
    : DefaultPrior(p)
{
    m_log = rundata.GetLogger();
    m_filename = p.options.find("image")->second;
    m_image = rundata.GetVoxelData(m_filename).AsRow();
}

void ImagePrior::DumpInfo(std::ostream &out) const
{
    out << "ImagePrior: Parameter " << m_idx << " '" << m_param_name << "'"
        << " filename: " << m_filename << " precision: " << m_params.prec();
}

double ImagePrior::ApplyToMVN(MVNDist *prior, const RunContext &ctx)
{
    prior->means(m_idx + 1) = m_image(ctx.v);

    SymmetricMatrix prec = prior->GetPrecisions();
    prec(m_idx + 1, m_idx + 1) = m_params.prec();
    prior->SetPrecisions(prec);

    return 0;
}

void ARDPrior::DumpInfo(std::ostream &out) const
{
    out << "ARDPrior: Parameter " << m_idx << " '" << m_param_name << "'"
        << " initial mean: " << m_params.mean() << " initial precision: " << m_params.prec();
}

double ARDPrior::ApplyToMVN(MVNDist *prior, const RunContext &ctx)
{
    SymmetricMatrix cov = prior->GetCovariance();
    double post_mean = ctx.fwd_post[ctx.v - 1].means(m_idx + 1);
    double post_cov = ctx.fwd_post[ctx.v - 1].GetCovariance()(m_idx + 1, m_idx + 1);
    // (Chappel et al 2009 Eq D4)
    double new_cov = post_mean * post_mean + post_cov;

    if (ctx.it == 0)
    {
        // Special case for first iter
        // Initially set prior to model default. The other alternative is to set
        // it to be initially non-informative, however this can still be achieved
        // by specifying the model prior mean/precision by options.
        cov(m_idx + 1, m_idx + 1) = m_params.var();
        prior->means(m_idx + 1) = m_params.mean();
        // LOG << "first iter ARD: " << m_params.var() << ", " << m_params.mean() << endl;
    }
    else
    {
        // LOG << "post: " << post_cov << ", " << post_mean << endl;
        // Update covariance on subsequent iterations
        cov(m_idx + 1, m_idx + 1) = new_cov;
        // LOG << "subs iter ARD: " << new_cov << ", " << prior->means(m_idx+1) << endl;
    }
    prior->SetCovariance(cov);

    // Calculate the free energy contribution from ARD term
    // (Chappel et al 2009, end of Appendix D)
    double b = 2 / new_cov;
    return -1.5 * (log(b) + digamma(0.5)) - 0.5 - gammaln(0.5) - 0.5 * log(b);
}

SpatialPrior::SpatialPrior(const Parameter &p, FabberRunData &rundata)
    : DefaultPrior(p)
    , m_aK(1e-8)
    , m_spatial_dims(3)
    , m_spatial_speed(-1)
{
    m_log = rundata.GetLogger();
    m_spatial_dims = rundata.GetIntDefault("spatial-dims", 3);
    if (m_spatial_dims < 0 || m_spatial_dims > 3)
    {
        throw InvalidOptionValue("spatial-dims", stringify(m_spatial_dims), "Must be 0, 1, 2 or 3");
    }
    else if (m_spatial_dims == 1)
    {
        WARN_ONCE("spatial-dims=1 is very weird... hope you're just testing!");
    }
    else if (m_spatial_dims == 2)
    {
        WARN_ONCE("spatial-dims=2 doesn't decompose into slices");
    }

    // FIXME check valid range
    m_spatial_speed = rundata.GetDoubleDefault("spatial-speed", -1);

    // FIXME still needed?
    m_update_first_iter = rundata.GetBool("update-spatial-prior-on-first-iteration");

    m_q1 = rundata.GetDoubleDefault("spatial-q1", 10.0);
    m_q2 = rundata.GetDoubleDefault("spatial-q2", 1.0);
}

void SpatialPrior::DumpInfo(std::ostream &out) const
{
    out << "SpatialPrior: Parameter " << m_idx << " '" << m_param_name << "'"
        << " type " << m_type_code << " mean: " << m_params.mean()
        << " precision: " << m_params.prec();
}

double SpatialPrior::CalculateaK(const RunContext &ctx)
{
    // Calculation of update equations in Penny et al 2005 Fig 4 (Spatial Precisions)
    //
    // sigmaK = Voxelwise parameter posterior variance
    // wK     = Voxelwise parameter posterior means
    // S      = Spatial precisions matrix. For Penny priors this is the Laplacian
    //          whereas for MRF priors we work with a matrix representing S' * S directly
    //
    // Theory notes and references are from MSC and should not be considered
    // reliable! Read Penny 2005 for details

    double trace_term = 0.0;    // First term for gk:   Tr(sigmaK*S'*S)
    double term2 = 0.0;         // Second term for gk:  wK'S'SwK
    for (int v = 1; v <= ctx.nvoxels; v++)
    {
        // Ignore voxels where numerical issues have occurred. Note that
        // excluded voxels are also deleted from neighbour lists for other voxels
        if (std::find(ctx.ignore_voxels.begin(), ctx.ignore_voxels.end(), v)
            != ctx.ignore_voxels.end()) continue;

        // Parameter variance
        double sigmaK = ctx.fwd_post.at(v - 1).GetCovariance()(m_idx + 1, m_idx + 1);

        // Number of neighbours
        int nn = ctx.neighbours.at(v - 1).size();

        if (m_type_code == PRIOR_SPATIAL_m)
        {
            // Markov random field without boundary correction
            // Assuming spatial_dims*2 nearest neighbours
            trace_term += sigmaK * m_spatial_dims * 2;
        }
        else if (m_type_code == PRIOR_SPATIAL_M)
        {
            // Markov random field with boundary correction
            // Using the actual number of nearest neighbours
            // (1e-8 term is to guarantee invertibility?)
            trace_term += sigmaK * (nn + 1e-8);
        }
        else if (m_type_code == PRIOR_SPATIAL_p)
        {
            // Penny prior without boundary correction
            // Uses Laplacian spatial matrix with
            // number of nearest neighbours = 2*spatial_dims
            trace_term += sigmaK * (4 * m_spatial_dims * m_spatial_dims + 2 * m_spatial_dims);
        }
        else
        {
            // Penny prior with boundary correction using actual
            // number of nearest neighbours
            trace_term += sigmaK * (nn * nn + nn);
        }

        // Posterior means
        double wK = ctx.fwd_post.at(v - 1).means(m_idx + 1);

        // Contribution from nearest neighbours - sum of differences
        // between voxel mean and neighbour mean
        double SwK = 0.0;
        for (vector<int>::const_iterator v2It = ctx.neighbours[v - 1].begin();
             v2It != ctx.neighbours.at(v - 1).end(); ++v2It)
        {
            SwK += wK - ctx.fwd_post.at(*v2It - 1).means(m_idx + 1);
        }

        // For priors with no boundary correction assume fixed number of neighbours
        // which means at boundaries some of the wK contribution will not have
        // appeared in the above sum. This is equivalent to assuming a mean
        // outside the boundary of zero (hence biased)
        if (m_type_code == PRIOR_SPATIAL_p || m_type_code == PRIOR_SPATIAL_m)
            SwK += wK * (m_spatial_dims * 2 - ctx.neighbours.at(v - 1).size());

        // For MRF spatial prior the spatial precision matrix S'S is handled
        // directly so we are effectively calculating wK * D * wK where
        // D is the spatial matrix. For Penny prior we work with the
        // Laplacian matrix and need wK * S' * S * wK
        if (m_type_code == PRIOR_SPATIAL_m || m_type_code == PRIOR_SPATIAL_M)
            term2 += SwK * wK;
        else
            term2 += SwK * SwK;
    }

    LOG << "SpatialPrior::Calculate aK " << m_idx << ": trace_term=" << trace_term << ", term2=" << term2 << endl;

    // Fig 4 in Penny (2005) update equations for gK, hK and aK
    //
    // Following Penny, prior on aK is a relatively uninformative gamma distribution with
    // q1 (theta) = 10 (1/q1 = 0.1) and q2 (k) = 1.0
    double gk = 1 / (0.5 * trace_term + 0.5 * term2 + 1/m_q1);
    double hK = (ctx.nvoxels * 0.5 + m_q2);
    double aK = gk * hK;

    if (aK < 1e-50)
    {
        // Don't let aK get too small
        LOG << "SpatialPrior::Calculate aK " << m_idx << ": was " << aK << endl;
        WARN_ONCE("SpatialPrior::Calculate aK - value was tiny - fixing to 1e-50");
        aK = 1e-50;
    }

    // Controls the speed of changes to aK - unsure whether this is useful or not but
    // it is only used if m_spatial_speed is given as an option

    // FIXME default for m_spatial_speed is -1 so code below will always be executed -
    // harmless but potentially confusing
    double aKMax = aK * m_spatial_speed;
    if (aKMax < 0.5)
    {
        // totally the wrong scaling.. oh well
        aKMax = 0.5;
    }

    if ((m_spatial_speed > 0) && (aK > aKMax))
    {
        LOG << "SpatialPrior::Calculate aK " << m_idx
            << ": Rate-limiting the increase on aK: was " << aK << ", now " << aKMax
            << endl;
        aK = aKMax;
    }

    LOG << "SpatialPrior::Calculate aK " << m_idx << ": New aK: " << aK << endl;
    return aK;
}

double SpatialPrior::ApplyToMVN(MVNDist *prior, const RunContext &ctx)
{
    // Comments and theory notes are from MSC and should not be trusted

    if (ctx.v == 1 && (ctx.it > 0 || m_update_first_iter))
    {
        // m_aK is a global (all voxels) spatial precision variable for
        // the parameter. It determines the degree of smoothness
        // and is estimated from the data. We update it on the first
        // voxel (unless it is the first iteration and m_update_first_iter
        // is false). See Penny et all 2004
        m_aK = CalculateaK(ctx);
    }

    // Loop over nearest neighbours of the current voxel
    // These have weighting +8
    int nn = ctx.neighbours[ctx.v - 1].size();
    double contrib_nn = 0.0;
    for (vector<int>::const_iterator nidIt = ctx.neighbours[ctx.v - 1].begin();
         nidIt != ctx.neighbours[ctx.v - 1].end(); ++nidIt)
    {
        int nid = *nidIt;
        const MVNDist &neighbourPost = ctx.fwd_post[nid - 1];
        contrib_nn += neighbourPost.means(m_idx + 1);
    }

    // Loop over second neighbours of the current voxel. Note that this list
    // contains duplicates of voxels which can be reached by two different routes.
    // By giving each a weighting of -1, this automatically generates the weightings
    // of -1 and -2 for linearly/diagonally connected voxels, as shown in Penny et
    // al 2004, Fig 3
    int nn2 = ctx.neighbours2[ctx.v - 1].size();
    double contrib_nn2 = 0.0;
    for (vector<int>::const_iterator nidIt = ctx.neighbours2[ctx.v - 1].begin();
         nidIt != ctx.neighbours2[ctx.v - 1].end(); ++nidIt)
    {
        int nid = *nidIt;
        const MVNDist &neighbourPost = ctx.fwd_post[nid - 1];
        contrib_nn2 += -neighbourPost.means(m_idx + 1);
    }

    // In priors without boundary correction, the number of neighbours is fixed by
    // the spatial dimensions, although this is not correct at the edges of the volume
    // The contributions from first and second neighbours are unchanged - this is
    // equivalent to the 'missing' voxels outside the boundary having value 0
    // (hence biased towards zero as noted in Penny 2005).
    if ((m_type_code == PRIOR_SPATIAL_p) || (m_type_code == PRIOR_SPATIAL_m))
    {
        assert(nn <= m_spatial_dims * 2);
        nn = 2 * m_spatial_dims;
        nn2 = 4 * m_spatial_dims * m_spatial_dims - nn;
    }

    // Prior precision
    //
    // The precision of the prior depends on the degree of spatial regularization
    // (m_aK) and the number of neighbours.

    // Calculate spatial precision for each prior type. Note the 1e-8 contribution
    // for the MRF prior, this is to ensure invertibility?
    double spatial_prec = 0;
    if (m_type_code == PRIOR_SPATIAL_M)
        spatial_prec = m_aK * (nn + 1e-8);
    else if (m_type_code == PRIOR_SPATIAL_m)
        spatial_prec = m_aK * nn;
    else if ((m_type_code == PRIOR_SPATIAL_P) || (m_type_code == PRIOR_SPATIAL_p))
        spatial_prec = m_aK * (nn * nn + nn);
    else
        assert(false);

    SymmetricMatrix precs = prior->GetPrecisions();
    if ((m_type_code == PRIOR_SPATIAL_p) || (m_type_code == PRIOR_SPATIAL_m))
    {
        // Penny-style DirichletBC prior ignores original prior precision completely
        precs(m_idx + 1, m_idx + 1) = spatial_prec;
    }
    else
    {
        // All other priors set a prior precision based on the sum of the original
        // prior precision and the spatial precision - i.e. spatial smoothing
        // makes prior more informative
        precs(m_idx + 1, m_idx + 1) = m_params.prec() + spatial_prec;
    }
    prior->SetPrecisions(precs);

    // Prior mean
    //
    // The prior mean depends on the degree of spatial regularization (which has been
    // absorbed into the prior precision/covariance matrix). It may also reflect the
    // original prior mean for the parameter, i.e. spatial regularization does not completely
    // override the original prior mean.

    // NB that we multiply by reciprocals rather than dividing. This is
    // to maximise numerical compatibility with NEWMAT which presumably
    // does it as an optimization when dividing a whole matrix by a constant
    double spatial_mean;
    if (m_type_code == PRIOR_SPATIAL_m)
    {
        // Dirichlet BCs on MRF i.e. no boundary correction
        double rec = 1 / double(nn);
        spatial_mean = contrib_nn * rec;
    }
    else if (m_type_code == PRIOR_SPATIAL_M)
    {
        double rec = 1 / double(nn);
        spatial_mean = contrib_nn * rec;
    }
    else if (nn != 0)
    {
        double rec = 1 / (8*nn - nn2);
        spatial_mean = (8*contrib_nn + contrib_nn2) * rec;
    }
    else
        spatial_mean = 0;

    if (m_type_code == PRIOR_SPATIAL_m || m_type_code == PRIOR_SPATIAL_M)
    {
        // These priors do not take account of the original parameter prior mean at all -
        // the prior mean is determined entirely by the expectation of spatial uniformity
        // compared to neighbouring voxels.
        //
        // For prior type m this reduces to spatial_mean since the covariance is in this
        // case simply the reciprocal of the spatial precision
        prior->means(m_idx + 1) = prior->GetCovariance()(m_idx + 1, m_idx + 1) * spatial_prec
            * spatial_mean;
    }
    else
    {
        // These priors include the original parameter prior mean as well as the values of
        // the parameter at neighbouring voxels.
        //
        // For prior type p, this reduces to spatial_mean + spatial_cov * original_prec * original_mean
        // For a weak original non-spatial prior (i.e. m_params->prec is small)  this reduces to
        // the spatial spatial_mean,
        prior->means(m_idx + 1) = prior->GetCovariance()(m_idx + 1, m_idx + 1)
            * (spatial_prec * spatial_mean + m_params.prec() * m_params.mean());
    }

    // Spatial prior does not contribute to the free energy?
    // Technically I think it should via the prior on ak (q1, q2) however this is quite uninformative.
    // It may not matter as we are not maximising F directly but rather following the update equations.
    return 0;
}

std::vector<Prior *> PriorFactory::CreatePriors(const std::vector<Parameter> &params)
{
    vector<Prior *> priors;
    for (size_t i = 0; i < params.size(); i++)
    {
        Prior *prior = PriorFactory::CreatePrior(params[i]);
        LOG << "PriorFactory::CreatePriors " << *prior << endl;
        priors.push_back(prior);
    }

    return priors;
}

PriorFactory::PriorFactory(FabberRunData &rundata)
    : Loggable(rundata.GetLogger())
    , m_rundata(rundata)
{
}

Prior *PriorFactory::CreatePrior(Parameter p)
{
    switch (p.prior_type)
    {
    case PRIOR_NORMAL:
    case PRIOR_DEFAULT:
        return new DefaultPrior(p);
    case PRIOR_IMAGE:
        return new ImagePrior(p, m_rundata);
    case PRIOR_SPATIAL_M:
    case PRIOR_SPATIAL_m:
    case PRIOR_SPATIAL_P:
    case PRIOR_SPATIAL_p:
        return new SpatialPrior(p, m_rundata);
    case PRIOR_ARD:
        return new ARDPrior(p, m_rundata);
    default:
        throw InvalidOptionValue("Prior type", stringify(p.prior_type), "Supported types: NMmPpAI");
    }
}
