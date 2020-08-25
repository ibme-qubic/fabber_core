/*  inference_vb.cc - Variational Bayes with optional spatial smoothing

 Adrian Groves and Matthew Webster, FMRIB Image Analysis Group

 Copyright (C) 2007-2010 University of Oxford  */

/*  CCOPYRIGHT */

#include "inference_vb.h"

#include "convergence.h"
#include "easylog.h"
#include "priors.h"
#include "run_context.h"
#include "tools.h"
#include "version.h"

#include <miscmaths/miscmaths.h>
#include <newmatio.h>

#include <math.h>

using MISCMATHS::sign;

static OptionSpec OPTIONS[] = {
    { "noise", OPT_STR, "Noise model to use (white or ar1)", OPT_REQ, "" },
    { "convergence", OPT_STR, "Name of method for detecting convergence", OPT_NONREQ, "maxits" },
    { "max-iterations", OPT_STR,
        "number of iterations of VB to use with the maxits convergence detector", OPT_NONREQ,
        "10" },
    { "min-fchange", OPT_STR,
        "When using the fchange convergence detector, the change in F to stop at", OPT_NONREQ,
        "10" },
    { "max-trials", OPT_STR, "When using the trial mode convergence detector, the maximum number "
                             "of trials after an initial reduction in F",
        OPT_NONREQ, "10" },
    { "print-free-energy", OPT_BOOL, "Output the free energy", OPT_NONREQ, "" },
    { "mcsteps", OPT_INT, "Number of motion correction steps", OPT_NONREQ, "0" },
    { "continue-from-mvn", OPT_MVN, "Continue previous run from output MVN files", OPT_NONREQ, "" },
    { "output-only", OPT_BOOL, "Skip model fitting, just output requested data based on supplied "
                               "MVN. Can only be used with continue-from-mvn",
        OPT_NONREQ, "" },
    { "noise-initial-prior", OPT_MATRIX, "MVN of initial noise prior", OPT_NONREQ, "" },
    { "noise-initial-posterior", OPT_MATRIX, "MVN of initial noise posterior", OPT_NONREQ, "" },
    { "noise-pattern", OPT_STR, "repeating pattern of noise variances for each point (e.g. 12 "
                                "gives odd and even data points different variances)",
        OPT_NONREQ, "1" },
    { "PSP_byname<n>", OPT_STR, "Name of model parameter to use image prior", OPT_NONREQ, "" },
    { "PSP_byname<n>_type", OPT_STR, "Type of image prior to use for parameter <n> - I=image prior",
        OPT_NONREQ, "" },
    { "PSP_byname<n>_image", OPT_IMAGE, "Image prior for parameter <n>", OPT_NONREQ, "" },
    { "PSP_byname<n>_prec", OPT_FLOAT, "Precision to apply to image prior for parameter <n>",
        OPT_NONREQ, "" },
    { "PSP_byname<n>_transform", OPT_STR, "Transform to apply to parameter <n>", OPT_NONREQ, "" },
    { "allow-bad-voxels", OPT_BOOL,
        "Continue if numerical error found in a voxel, rather than stopping", OPT_NONREQ, "" },
    { "ar1-cross-terms", OPT_STR, "For AR1 noise, type of cross-linking (dual, same or none)",
        OPT_NONREQ, "dual" },
    { "spatial-dims", OPT_INT, "Number of spatial dimensions", OPT_NONREQ, "3" },
    { "spatial-speed", OPT_STR, "Restrict speed of spatial smoothing", OPT_NONREQ, "-1" },
    { "distance-measure", OPT_STR, "", OPT_NONREQ, "dist1" },
    { "param-spatial-priors", OPT_STR,
        "Type of spatial priors for each parameter, as a sequence of characters. "
        "N=nonspatial, M=Markov random field, P=Penny, A=ARD",
        OPT_NONREQ, "N+" },
    { "update-spatial-prior-on-first-iteration", OPT_BOOL, "", OPT_NONREQ, "" },
    { "locked-linear-from-mvn", OPT_MVN, "MVN file containing fixed centres for linearization",
        OPT_NONREQ, "" },
    { "" },
};

void Vb::GetOptions(vector<OptionSpec> &opts) const
{
    InferenceTechnique::GetOptions(opts);

    for (int i = 0; OPTIONS[i].name != ""; i++)
    {
        opts.push_back(OPTIONS[i]);
    }
}

string Vb::GetDescription() const
{
    return "Variational Bayes inference technique";
}
string Vb::GetVersion() const
{
    return fabber_version();
}
InferenceTechnique *Vb::NewInstance()
{
    return new Vb();
}
void Vb::Initialize(FwdModel *fwd_model, FabberRunData &rundata)
{
    InferenceTechnique::Initialize(fwd_model, rundata);

    // Get noise model.
    m_noise = std::auto_ptr<NoiseModel>(NoiseModel::NewFromName(rundata.GetString("noise")));
    m_noise->Initialize(rundata);
    m_noise_params = m_noise->NumParams();
    LOG << "Vb::Noise has " << m_noise_params << " parameters" << endl;

    // Figure out if F needs to be calculated every iteration
    m_saveF = rundata.GetBool("save-free-energy");
    m_saveFsHistory = rundata.GetBool("save-free-energy-history");
    m_printF = rundata.GetBool("print-free-energy");

    // Motion correction related setup - by default no motion correction
    m_num_mcsteps = convertTo<int>(rundata.GetStringDefault("mcsteps", "0"));

    m_spatial_dims = rundata.GetIntDefault("spatial-dims", 3, 0, 3);
    if (m_spatial_dims == 1)
    {
        WARN_ONCE("spatial-dims=1 is weird... hope you're just testing!");
    }
    else if (m_spatial_dims == 2)
    {
        WARN_ONCE("spatial-dims=2 may not work the way you expect");
    }

    // Locked linearizations, if requested
    m_locked_linear = rundata.GetStringDefault("locked-linear-from-mvn", "") != "";
}

void Vb::InitializeNoiseFromParam(FabberRunData &rundata, NoiseParams *dist, string param_key)
{
    string filename = rundata.GetStringDefault(param_key, "modeldefault");
    if (filename != "modeldefault")
    {
        // FIXME should there be checking of size here
        LOG << "VbInferenceTechnique::Loading " << param_key << " distribution from " << filename
            << endl;
        dist->InputFromMVN(MVNDist(filename, m_log));
    }
}

void Vb::SetupPerVoxelDists(FabberRunData &rundata)
{
    // Initialized in voxel loop below (from file or default as required)
    m_ctx->noise_post.resize(m_nvoxels, NULL);
    m_ctx->noise_prior.resize(m_nvoxels, NULL);
    m_ctx->fwd_post.resize(m_nvoxels);

    // Re-centred in voxel loop below
    m_lin_model.resize(m_nvoxels, LinearizedFwdModel(m_model));

    // Initialized in voxel loop below
    m_conv.resize(m_nvoxels, NULL);
    string conv_name = rundata.GetStringDefault("convergence", "maxits");

    // Model prior is updated during main voxel loop
    m_ctx->fwd_prior.resize(m_nvoxels, MVNDist(m_num_params, m_log));

    // Loaded from file if required, otherwise initialized during calculation
    resultMVNs.resize(m_nvoxels, NULL);

    // Initialized during calculation
    resultFs.resize(m_nvoxels, 9999); // 9999 is a garbage default value
    resultFsHistory.resize(m_nvoxels);

    // Whether to fix the linearization centres (default: false)
    vector<MVNDist *> lockedLinearDists;
    Matrix lockedLinearCentres;
    if (m_locked_linear)
    {
        string file = rundata.GetString("locked-linear-from-mvn");
        LOG << "Vb::Loading fixed linearization centres from the MVN '" << file
            << "'\nNOTE: This does not check if the correct "
               "number of parameters is present!\n";
        MVNDist::Load(lockedLinearDists, file, rundata, m_log);
        lockedLinearCentres.ReSize(m_num_params, m_nvoxels);
    }

    // If we are resuming from a previous run, there will be data containing a per-voxel
    // distribution of the model parameters, and noise as well.
    bool continueFromMvn = false;
    try {
        rundata.GetVoxelData("continue-from-mvn");
        continueFromMvn = true;
    }
    catch(DataNotFound &e) {
        // no worries
    }

    if (continueFromMvn)
    {
        LOG << "Vb::Continuing from MVN" << endl;
        // Optional list of parameters in MVN
        string paramFilename = rundata.GetStringDefault("continue-from-params", "");
        InitMVNFromFile(rundata, paramFilename);
    }

    // Initial noise distributions
    auto_ptr<NoiseParams> initialNoisePrior(m_noise->NewParams());
    auto_ptr<NoiseParams> initialNoisePosterior(m_noise->NewParams());
    m_noise->HardcodedInitialDists(*initialNoisePrior, *initialNoisePosterior);
    InitializeNoiseFromParam(rundata, initialNoisePrior.get(), "noise-initial-prior");
    InitializeNoiseFromParam(rundata, initialNoisePosterior.get(), "noise-initial-posterior");

    for (int v = 1; v <= m_nvoxels; v++)
    {
        if (continueFromMvn)
        {
            m_ctx->fwd_post[v - 1] = resultMVNs.at(v - 1)->GetSubmatrix(1, m_num_params);
            assert(m_num_params + m_noise_params == resultMVNs.at(v - 1)->GetSize());
            m_ctx->noise_post[v - 1] = m_noise->NewParams();
            m_ctx->noise_post[v - 1]->InputFromMVN(resultMVNs.at(v - 1)->GetSubmatrix(
                m_num_params + 1, m_num_params + m_noise_params));
        }
        else
        {
            // Set the initial posterior for model params. Model
            // may want the voxel data in order to do this
            PassModelData(v);
            m_model->GetInitialPosterior(m_ctx->fwd_post[v - 1], rundata);
            // Set initial noise posterior
            m_ctx->noise_post[v - 1] = initialNoisePosterior->Clone();
        }

        if (m_locked_linear)
        {
            lockedLinearCentres.Column(v)
                = lockedLinearDists.at(v - 1)->means.Rows(1, m_num_params);
            m_lin_model[v - 1].ReCentre(lockedLinearCentres.Column(v));
        }
        else
        {
            m_lin_model[v - 1].ReCentre(m_ctx->fwd_post[v - 1].means);
        }

        // Create per-voxel convergence detector. Initialization of m_needF is
        // inefficient but not harmful because all convergence detectors are the same type
        m_conv[v - 1] = ConvergenceDetector::NewFromName(conv_name);
        m_conv[v - 1]->Initialize(rundata);
        m_needF = m_conv[v - 1]->UseF() || m_printF || m_saveF || m_saveFsHistory;

        m_ctx->noise_prior[v - 1] = initialNoisePrior->Clone();
        m_noise->Precalculate(
            *m_ctx->noise_post[v - 1], *m_ctx->noise_prior[v - 1], m_origdata->Column(v));
    }
}

void Vb::PassModelData(int v)
{
    // Pass in data, coords and supplemental data for this voxel
    ColumnVector data = m_origdata->Column(v);
    ColumnVector vcoords = m_coords->Column(v);
    if (m_suppdata->Ncols() > 0)
    {
        ColumnVector suppy = m_suppdata->Column(v);
        m_model->PassData(v, data, vcoords, suppy);
    }
    else
    {
        m_model->PassData(v, data, vcoords);
    }
}

void Vb::IgnoreVoxel(int v)
{
    LOG << "Vb::IgnoreVoxel This voxel will be ignored in further updates" << endl;

    m_ctx->ignore_voxels.push_back(v);

    // Remove voxel from lists of neighbours of other voxels.
    // We identify affected voxels by looking in the neighbour
    // lists for the bad voxel, because any voxel which has
    // the bad voxel as a neighbour will be a neighbour of the
    // bad voxel
    vector<int> nn = m_ctx->neighbours[v - 1];
    for (vector<int>::iterator i = nn.begin(); i != nn.end(); ++i)
    {
        // Reference to list of neighbours of some other voxel which
        // has the bad voxel as a neighbour
        vector<int> &n2 = m_ctx->neighbours[*i - 1];

        n2.erase(std::remove(n2.begin(), n2.end(), v), n2.end());
    }

    // Same for next-nearest-neighbours
    nn = m_ctx->neighbours2[v - 1];
    for (vector<int>::iterator i = nn.begin(); i != nn.end(); ++i)
    {
        // Reference to list of neighbours of some other voxel which
        // has the bad voxel as a neighbour
        vector<int> &n2 = m_ctx->neighbours2[*i - 1];

        n2.erase(std::remove(n2.begin(), n2.end(), v), n2.end());
    }
}

/**
 * Calculate free energy. Note that this is currently unused in spatial VB
 */
double Vb::CalculateF(int v, string label, double Fprior)
{
    double F = 1234.5678;
    if (m_needF)
    {
        F = m_noise->CalcFreeEnergy(*m_ctx->noise_post[v - 1], *m_ctx->noise_prior[v - 1],
            m_ctx->fwd_post[v - 1], m_ctx->fwd_prior[v - 1], m_lin_model[v - 1],
            m_origdata->Column(v));
        F += Fprior;
        resultFs[v - 1] = F;
        if (m_printF)
        {
            LOG << "Vb::F" << label << " = " << F << endl;
        }
    }
    return F;
}

void Vb::DebugVoxel(int v, const string &where)
{
    LOG << where << " - voxel " << v << " of " << m_nvoxels << endl;
    LOG << "Prior means: " << endl << m_ctx->fwd_prior[v - 1].means.t();
    LOG << "Prior precisions: " << endl << m_ctx->fwd_prior[v - 1].GetPrecisions();
    LOG << "Posterior means: " << endl << m_ctx->fwd_post[v - 1].means.t();
    LOG << "Noise prior means: " << endl << m_ctx->noise_prior[v - 1]->OutputAsMVN().means.t();
    LOG << "Noise prior precisions: " << endl
        << m_ctx->noise_prior[v - 1]->OutputAsMVN().GetPrecisions();
    LOG << "Centre: " << endl << m_lin_model[v - 1].Centre();
    LOG << "Offset: " << endl << m_lin_model[v - 1].Offset();
    LOG << "Jacobian: " << endl << m_lin_model[v - 1].Jacobian() << endl;
}

bool Vb::IsSpatial(FabberRunData &rundata) const
{
    if (rundata.GetString("method") == "spatialvb")
    {
        return true;
    }
    else
    {
        // Really clunky way to detect if any spatial priors have been specified
        vector<Parameter> params;
        m_model->GetParameters(rundata, params);
        for (vector<Parameter>::iterator iter = params.begin(); iter != params.end(); iter++)
        {
            switch (iter->prior_type)
            {
            case PRIOR_SPATIAL_M:
            case PRIOR_SPATIAL_m:
            case PRIOR_SPATIAL_P:
            case PRIOR_SPATIAL_p:
                return true;
            }
        }
    }
    return false;
}

void Vb::DoCalculations(FabberRunData &rundata)
{
    // extract data (and the coords) from rundata for the (first) VB run
    // Rows are volumes
    // Columns are (time) series
    // num Rows is size of (time) series
    // num Cols is size of volumes
    m_origdata = &rundata.GetMainVoxelData();
    m_coords = &rundata.GetVoxelCoords();
    m_suppdata = &rundata.GetVoxelSuppData();
    m_nvoxels = m_origdata->Ncols();
    m_ctx = new RunContext(m_nvoxels);

    // pass in some (dummy) data/coords here just in case the model relies upon it
    // use the first voxel values as our dummies FIXME this shouldn't really be
    // necessary, need to find way for model to know about the data beforehand.
    if (m_nvoxels > 0)
        PassModelData(1);

    // Only call DoCalculations once
    assert(resultMVNs.empty());
    assert(resultFs.empty());

    SetupPerVoxelDists(rundata);

    if (rundata.GetBool("output-only"))
    {
        // Do no calculations - now we have set resultMVNs we can finish
        LOG << "Vb::DoCalculations output-only set - not performing any calculations" << endl;
    }
    else if (IsSpatial(rundata))
    {
        DoCalculationsSpatial(rundata);
    }
    else
    {
        DoCalculationsVoxelwise(rundata);
    }

    if (!m_needF)
    {
        // clearing resultFs here should prevent an F image from being saved.
        resultFs.clear();
    }

    // Delete stuff (avoid memory leaks)
    for (int v = 1; v <= m_nvoxels; v++)
    {
        delete m_ctx->noise_post[v - 1];
        delete m_ctx->noise_prior[v - 1];
        delete m_conv[v - 1];
    }
    delete m_ctx;
}

void Vb::DoCalculationsVoxelwise(FabberRunData &rundata)
{
    vector<Parameter> params;
    m_model->GetParameters(rundata, params);
    vector<Prior *> priors = PriorFactory(rundata).CreatePriors(params);

    LOG << "Vb::Voxelwise calculations loop" << endl;
    // Loop over voxels
    for (int v = 1; v <= m_nvoxels; v++)
    {
        PassModelData(v);

        m_ctx->v = v;
        m_ctx->it = 0;

        // Save our model parameters in case we need to revert later.
        // Note need to save prior in case ARD is being used
        NoiseParams *const noisePosteriorSave = m_ctx->noise_post[v - 1]->Clone();
        MVNDist fwdPosteriorSave(m_ctx->fwd_post[v - 1]);
        MVNDist fwdPriorSave(m_ctx->fwd_prior[v - 1]);

        // Give an indication of the progress through the voxels;
        rundata.Progress(v, m_nvoxels);
        double F = 1234.5678;
        double Fprior = 0;

        try
        {
            m_lin_model[v - 1].ReCentre(m_ctx->fwd_post[v - 1].means);
            m_conv[v - 1]->Reset();

            // START the VB updates and run through the relevant iterations (according to the
            // convergence testing)
            do
            {
                // Save old values if the convergence detector found that they were the best so far
                if (m_conv[v - 1]->NeedSave())
                {
                    *noisePosteriorSave = *m_ctx->noise_post[v - 1]; // copy values, not pointer!
                    fwdPosteriorSave = m_ctx->fwd_post[v - 1];
                    fwdPriorSave = m_ctx->fwd_prior[v - 1];
                    if (m_debug)
                        DebugVoxel(v, "Saving as best solution so far");
                }

                for (int k = 0; k < m_num_params; k++)
                {
                    Fprior = priors[k]->ApplyToMVN(&m_ctx->fwd_prior[v - 1], *m_ctx);
                }

                if (m_debug)
                    DebugVoxel(v, "Applied priors");

                F = CalculateF(v, "before", Fprior);

                m_noise->UpdateTheta(*m_ctx->noise_post[v - 1], m_ctx->fwd_post[v - 1],
                    m_ctx->fwd_prior[v - 1], m_lin_model[v - 1], m_origdata->Column(v), NULL,
                    m_conv[v - 1]->LMalpha());

                if (m_debug)
                    DebugVoxel(v, "Updated params");

                F = CalculateF(v, "theta", Fprior);

                m_noise->UpdateNoise(*m_ctx->noise_post[v - 1], *m_ctx->noise_prior[v - 1],
                    m_ctx->fwd_post[v - 1], m_lin_model[v - 1], m_origdata->Column(v));

                if (m_debug)
                    DebugVoxel(v, "Updated noise");

                F = CalculateF(v, "phi", Fprior);

                // Linearization update
                // Update the linear model before doing Free energy calculation
                // (and ready for next round of theta and phi updates)
                m_lin_model[v - 1].ReCentre(m_ctx->fwd_post[v - 1].means);

                if (m_debug)
                    DebugVoxel(v, "Re-centered");

                F = CalculateF(v, "lin", Fprior);
                if (m_saveFsHistory)
                    resultFsHistory.at(v - 1).push_back(F);

                ++m_ctx->it;
            } while (!m_conv[v - 1]->Test(F));

            if (m_debug)
                LOG << "Converged after " << m_ctx->it << " iterations" << endl;

            // Save old values if best so far FIXME is this needed?
            if (m_conv[v - 1]->NeedSave())
            {
                *noisePosteriorSave = *m_ctx->noise_post[v - 1]; // copy values, not pointer!
                fwdPosteriorSave = m_ctx->fwd_post[v - 1];
                fwdPriorSave = m_ctx->fwd_prior[v - 1];
                if (m_debug)
                    DebugVoxel(v, "Saving as best solution at end");
            }

            // Revert to previous best values at last stage if required
            if (m_conv[v - 1]->NeedRevert())
            {
                *m_ctx->noise_post[v - 1] = *noisePosteriorSave;
                m_ctx->fwd_post[v - 1] = fwdPosteriorSave;
                m_ctx->fwd_prior[v - 1] = fwdPriorSave;
                m_lin_model[v - 1].ReCentre(m_ctx->fwd_post[v - 1].means);
                if (m_debug)
                        DebugVoxel(v, "Reverted to better solution");
                F = CalculateF(v, "revert", Fprior);
            }

            delete noisePosteriorSave;
        }
        catch (FabberInternalError &e)
        {
            LOG << "Vb::Internal error for voxel " << v << " at " << m_coords->Column(v).t()
                << " : " << e.what() << endl;

            if (m_halt_bad_voxel)
                throw;
        }
        catch (NEWMAT::Exception &e)
        {
            LOG << "Vb::NEWMAT exception for voxel " << v << " at " << m_coords->Column(v).t()
                << " : " << e.what() << endl;

            if (m_halt_bad_voxel)
                throw;
        }

        // now write the results to resultMVNs
        try
        {
            resultMVNs.at(v - 1)
                = new MVNDist(m_ctx->fwd_post[v - 1], m_ctx->noise_post[v - 1]->OutputAsMVN());
            if (m_needF)
                resultFs.at(v - 1) = F;
            if (m_saveFsHistory)
                resultFsHistory.at(v - 1).push_back(F);
        }
        catch (...)
        {
            // Even that can fail, due to results being singular
            LOG << "Vb::Can't give any sensible answer for this voxel; outputting zero +- "
                   "identity\n";
            MVNDist *tmp = new MVNDist(m_log);
            tmp->SetSize(m_ctx->fwd_post[v - 1].means.Nrows()
                + m_ctx->noise_post[v - 1]->OutputAsMVN().means.Nrows());
            tmp->SetCovariance(IdentityMatrix(tmp->means.Nrows()));
            resultMVNs.at(v - 1) = tmp;
            if (m_needF)
                resultFs.at(v - 1) = F;
            if (m_saveFsHistory)
                resultFsHistory.at(v - 1).push_back(F);
        }
    }
    for (unsigned int i = 0; i < priors.size(); i++)
    {
        delete priors[i];
    }
}

void Vb::DoCalculationsSpatial(FabberRunData &rundata)
{
    // Pass in some (dummy) data/coords here just in case the model relies upon it
    // use the first voxel values as our dummies FIXME this shouldn't really be
    // necessary, need to find way for model to know about the data beforehand.
    if (m_nvoxels > 0)
        PassModelData(1);

    // Make the neighbours[] lists if required
    // if (m_prior_types_str.find_first_of("mMpP") != string::npos)
    if (true) // FIXME
    {
        CalcNeighbours(*m_coords);
    }

    vector<Parameter> params;
    m_model->GetParameters(rundata, params);
    vector<Prior *> priors = PriorFactory(rundata).CreatePriors(params);

    // Spatial loop currently uses a global convergence detector FIXME
    // needs to change
    CountingConvergenceDetector conv;
    conv.Initialize(rundata);
    double Fglobal = 1234.5678;
    int maxits = convertTo<int>(rundata.GetStringDefault("max-iterations", "10"));

    // MAIN ITERATION LOOP
    do
    {
        LOG << endl << "*** Spatial iteration *** " << (m_ctx->it + 1) << endl;

        // Give an indication of the progress through the voxels;
        rundata.Progress(m_ctx->it, maxits);
        double Fprior = 0;

        // ITERATE OVER VOXELS
        for (int v = 1; v <= m_nvoxels; v++)
        {
            m_ctx->v = v;

            PassModelData(v);

            // The steps below are essentially the same as regular VB, although
            // the code looks different as the per-voxel dists are set up at the
            // start rather than as we go
            try
            {
                Fprior = 0;

                // Apply prior updates for spatial or ARD priors
                for (int k = 0; k < m_num_params; k++)
                {
                    Fprior += priors[k]->ApplyToMVN(&m_ctx->fwd_prior[v - 1], *m_ctx);
                }
                if (m_debug)
                    DebugVoxel(v, "Priors set");

                // Ignore voxels where numerical issues have occurred
                if (std::find(m_ctx->ignore_voxels.begin(), m_ctx->ignore_voxels.end(), v)
                    != m_ctx->ignore_voxels.end())
                {
                    LOG << "Ignoring voxel " << v << endl;
                    continue;
                }

                CalculateF(v, "before", Fprior);

                m_noise->UpdateTheta(*m_ctx->noise_post[v - 1], m_ctx->fwd_post[v - 1],
                    m_ctx->fwd_prior[v - 1], m_lin_model[v - 1], m_origdata->Column(v), NULL, 0);
                if (m_debug)
                    DebugVoxel(v, "Theta updated");

                CalculateF(v, "theta", Fprior);
            }
            catch (FabberInternalError &e)
            {
                LOG << "Vb::Internal error for voxel " << v << " at " << m_coords->Column(v).t()
                    << " : " << e.what() << endl;

                if (m_halt_bad_voxel)
                    throw;
                else
                    IgnoreVoxel(v);
            }
            catch (NEWMAT::Exception &e)
            {
                LOG << "Vb::NEWMAT exception for voxel " << v << " at " << m_coords->Column(v).t()
                    << " : " << e.what() << endl;

                if (m_halt_bad_voxel)
                    throw;
                else
                    IgnoreVoxel(v);
            }
        }

        Fglobal = 0;
        for (int v = 1; v <= m_nvoxels; v++)
        {
            try {
                // Ignore voxels where numerical issues have occurred
                if (std::find(m_ctx->ignore_voxels.begin(), m_ctx->ignore_voxels.end(), v)
                    != m_ctx->ignore_voxels.end())
                {
                    LOG << "Ignoring voxel " << v << endl;
                    continue;
                }

                PassModelData(v);

                m_noise->UpdateNoise(*m_ctx->noise_post[v - 1], *m_ctx->noise_prior[v - 1],
                    m_ctx->fwd_post[v - 1], m_lin_model[v - 1], m_origdata->Column(v));
                if (m_debug)
                    DebugVoxel(v, "Noise updated");

                CalculateF(v, "noise", Fprior);

                if (!m_locked_linear)
                    m_lin_model[v - 1].ReCentre(m_ctx->fwd_post[v - 1].means);
                if (m_debug)
                    DebugVoxel(v, "Re-centre");

                Fglobal += CalculateF(v, "lin", Fprior);
            }
            catch (FabberInternalError &e)
            {
                LOG << "Vb::Internal error for voxel " << v << " at " << m_coords->Column(v).t()
                    << " : " << e.what() << endl;

                if (m_halt_bad_voxel)
                    throw;
                else
                    IgnoreVoxel(v);
            }
            catch (NEWMAT::Exception &e)
            {
                LOG << "Vb::NEWMAT exception for voxel " << v << " at " << m_coords->Column(v).t()
                    << " : " << e.what() << endl;

                if (m_halt_bad_voxel)
                    throw;
                else
                    IgnoreVoxel(v);
            }
        }

        ++m_ctx->it;
    } while (!conv.Test(Fglobal));

    // Interesting addition: calculate "coefficient resels" from Penny et al. 2005
    for (int k = 1; k <= m_num_params; k++)
    {
        ColumnVector gamma_vk(m_nvoxels);
        for (int v = 1; v <= m_nvoxels; v++)
        {
            // Ignore voxels where numerical issues have occurred
            if (std::find(m_ctx->ignore_voxels.begin(), m_ctx->ignore_voxels.end(), v)
                != m_ctx->ignore_voxels.end())
            {
                gamma_vk(v) = 0;
                continue;
            }

            try
            {
                gamma_vk(v) = 1
                    - m_ctx->fwd_post[v - 1].GetCovariance()(k, k)
                        / m_ctx->fwd_prior[v - 1].GetCovariance()(k, k);
            }
            catch (...)
            {
                // Even that can fail, due to results being singular
                LOG << "Vb::Coefficient resels failed for voxel " << v << endl;
                gamma_vk(v) = 0;
            }
        }
        LOG << "Vb::Coefficient resels per voxel for param " << k << ": "
            << gamma_vk.Sum() / m_nvoxels << endl;
    }

    for (int v = 1; v <= m_nvoxels; v++)
    {
        resultMVNs[v - 1]
            = new MVNDist(m_ctx->fwd_post[v - 1], m_ctx->noise_post[v - 1]->OutputAsMVN());
    }
    for (unsigned int i = 0; i < priors.size(); i++)
    {
        delete priors[i];
    }
}

void Vb::CheckCoordMatrixCorrectlyOrdered(const Matrix &coords)
{
    // Only 3D
    assert(coords.Nrows() == 3);

    // Voxels are stored one per column, each column is the x/y/z coords
    const int m_nvoxels = coords.Ncols();

    // Go through each voxel one at a time apart from last
    for (int v = 1; v <= m_nvoxels - 1; v++)
    {
        // Find difference between current coords and next
        ColumnVector diff = coords.Column(v + 1) - coords.Column(v);

        // Check order
        // +1 = +x, +10 = +y, +100 = +z, -100 = -z+x, etc.
        int d = sign(diff(1)) + 10 * sign(diff(2)) + 100 * sign(diff(3));
        if (d <= 0)
        {
            LOG << "Vb::Found mis-ordered voxels " << v << " and " << v + 1 << ": d=" << d << endl;
            throw FabberInternalError(
                "Coordinate matrix must be in correct order to use adjacency-based priors.");
        }
    }
}

// Binary search for data(index) == num
// Assumes data is sorted ascending!!
// Either returns an index such that data(index) == num
//   or -1 if num is not present in data.
static inline int binarySearch(const ColumnVector &data, int num)
{
    int first = 1, last = data.Nrows();

    while (first <= last)
    {
        int test = (first + last) / 2;

        if (data(test) < num)
        {
            first = test + 1;
        }
        else if (data(test) > num)
        {
            last = test - 1;
        }
        else if (data(test) == num)
        {
            return test;
        }
        else
        {
            assert(false); // logic error!  data wasn't sorted?
        }
    }
    return -1;
}

/**
 * Calculate nearest and second-nearest neighbours for the voxels
 */
void Vb::CalcNeighbours(const Matrix &coords)
{
    const int nVoxels = coords.Ncols();
    if (nVoxels == 0)
        return;

    // Voxels must be ordered by increasing z, y and x values respectively
    // otherwise binary search for voxel by offset will not work
    CheckCoordMatrixCorrectlyOrdered(coords);

    // Create a column vector with one entry per voxel.
    ColumnVector offsets(nVoxels);

    // Populate offsets with the offset into the
    // matrix of each voxel. We assume that co-ordinates
    // could be zero but not negative
    int xsize = coords.Row(1).Maximum() + 1;
    int ysize = coords.Row(2).Maximum() + 1;
    for (int v = 1; v <= nVoxels; v++)
    {
        int x = coords(1, v);
        int y = coords(2, v);
        int z = coords(3, v);
        int offset = z * xsize * ysize + y * xsize + x;
        offsets(v) = offset;
    }

    // Delta is a list of offsets to find nearest
    // neighbours in x y and z direction (not diagonally)
    // Of course applying these offsets naively would not
    // always work, e.g. offset of -1 in the x direction
    // will not be a nearest neighbour for the first voxel
    // so need to check for this in subsequent code
    vector<int> delta;
    delta.push_back(1);              // next row
    delta.push_back(-1);             // prev row
    delta.push_back(xsize);          // next column
    delta.push_back(-xsize);         // prev column
    delta.push_back(xsize * ysize);  // next slice
    delta.push_back(-xsize * ysize); // prev slice

    // Don't look for neighbours in all dimensions.
    // For example if spatialDims=2, max_delta=3 so we
    // only look for neighbours in rows and columns
    //
    // However note we still need the full list of 3D deltas for later
    int max_delta = m_spatial_dims * 2 - 1;

    // Neighbours is a vector of vectors, so each voxel
    // will have an entry which is a vector of its neighbours
    m_ctx->neighbours.resize(nVoxels);

    // Go through each voxel. Note that offsets is indexed from 1 not 0
    // however the offsets themselves (potentially) start at 0.
    for (int vid = 1; vid <= nVoxels; vid++)
    {
        // Get the voxel offset into the matrix
        int pos = int(offsets(vid));

        // Now search for neighbours
        for (int n = 0; n <= max_delta; n++)
        {
            // is there a voxel at this neighbour position?
            // indexed from 1; id == -1 if not found.
            int id = binarySearch(offsets, pos + delta[n]);

            // No such voxel: continue
            if (id < 0)
                continue;

            // Check for wrap-around

            // Don't check for wrap around on final co-ord
            // PREVIOUSLY		if (delta.size() >= n + 2)
            // Changed (fixed)? because if spatialDims != 3 we still need
            // to check for wrap around in y-coordinate FIXME check
            if (n < 4)
            {
                bool ignore = false;
                if (delta[n] > 0)
                {
                    int test = delta[n + 2];
                    if (test > 0)
                        ignore = (pos % test) >= test - delta[n];
                }
                else
                {
                    int test = -delta[n + 2];
                    if (test > 0)
                        ignore = (pos % test) < -delta[n];
                }
                if (ignore)
                {
                    continue;
                }
            }

            // If we get this far, add it to the list
            m_ctx->neighbours.at(vid - 1).push_back(id);
        }
    }

    // Similar algorithm but looking for Neighbours-of-neighbours, excluding self,
    // but including duplicates if there are two routes to get there
    // (diagonally connected)
    m_ctx->neighbours2.resize(nVoxels);

    for (int vid = 1; vid <= nVoxels; vid++)
    {
        // Go through the list of neighbours for each voxel.
        for (unsigned n1 = 0; n1 < m_ctx->neighbours.at(vid - 1).size(); n1++)
        {
            // n1id is the voxel index (not the offset) of the neighbour
            int n1id = m_ctx->neighbours[vid - 1].at(n1);
            int checkNofN = 0;
            // Go through each of it's neighbours. Add each, apart from original voxel
            for (unsigned n2 = 0; n2 < m_ctx->neighbours.at(n1id - 1).size(); n2++)
            {
                int n2id = m_ctx->neighbours[n1id - 1].at(n2);
                if (n2id != vid)
                {
                    m_ctx->neighbours2[vid - 1].push_back(n2id);
                }
                else
                    checkNofN++;
            }

            if (checkNofN != 1)
            {
                throw FabberInternalError("Each of this voxel's neighbours must have "
                                          "this voxel as a neighbour");
            }
        }
    }
}

void Vb::SaveResults(FabberRunData &rundata) const
{
    InferenceTechnique::SaveResults(rundata);

    LOG << "Vb::Preparing to save results..." << endl;
    int nVoxels = resultMVNs.size();

    if (rundata.GetBool("save-noise-mean") | rundata.GetBool("save-noise-std"))
    {
        if (m_noise_params > 0)
        {
            LOG << "Vb::Writing noise" << endl;
            Matrix noiseMean, noiseStd;
            noiseMean.ReSize(m_noise_params, nVoxels);
            noiseStd.ReSize(m_noise_params, nVoxels);
            for (int vox = 1; vox <= nVoxels; vox++)
            {
                for (int i = 1; i <= m_noise_params; i++)
                {
                    noiseStd(i, vox) = sqrt(
                        resultMVNs[vox - 1]->GetCovariance()(i + m_num_params, i + m_num_params));
                    noiseMean(i, vox) = resultMVNs[vox - 1]->means(i + m_num_params);
                }
            }
            // FIXME was this being saved before? Should it be?
            if (rundata.GetBool("save-noise-mean"))
                rundata.SaveVoxelData("noise_means", noiseMean);
            if (rundata.GetBool("save-noise-std"))
                rundata.SaveVoxelData("noise_stdevs", noiseStd);
        }
    }

    // Save the Free Energy estimates
    if (m_saveF && !resultFs.empty())
    {
        LOG << "Vb::Writing free energy" << endl;
        assert((int)resultFs.size() == nVoxels);
        Matrix freeEnergy;
        freeEnergy.ReSize(1, nVoxels);
        for (int vox = 1; vox <= nVoxels; vox++)
        {
            freeEnergy(1, vox) = resultFs.at(vox - 1);
        }
        rundata.SaveVoxelData("freeEnergy", freeEnergy);
    }
    else
    {
        LOG << "Vb::Free energy wasn't recorded, so no freeEnergy data saved" << endl;
    }

    // Save the Free Energy history
    if ((nVoxels > 0)  && m_saveFsHistory && !resultFsHistory.empty())
    {
        LOG << "Vb::Writing free energy history" << endl;
        assert((int)resultFsHistory.size() == nVoxels);
        size_t num_iters = resultFsHistory[0].size();
        for (int vox = 1; vox <= nVoxels; vox++)
        {
            if (resultFsHistory[vox-1].size() > num_iters)
            {
                num_iters = resultFsHistory[vox-1].size();
            }
        }

        LOG << "Vb::Maximum number of iterations for a voxel was: " << num_iters << endl;
        Matrix freeEnergyHistory;
        freeEnergyHistory.ReSize(num_iters, nVoxels);
        for (int vox = 1; vox <= nVoxels; vox++)
        {
            for (unsigned int iter=1; iter <= num_iters; iter++)
            {
                if (iter <= resultFsHistory.at(vox - 1).size())
                {
                    freeEnergyHistory(iter, vox) = resultFsHistory.at(vox - 1).at(iter-1);
                }
                else
                {
                    freeEnergyHistory(iter, vox) = resultFsHistory.at(vox - 1).at(resultFsHistory.at(vox - 1).size()-1);
                }
            }
        }
        rundata.SaveVoxelData("freeEnergyHistory", freeEnergyHistory);
    }

    LOG << "Vb::Done writing results." << endl;
}
