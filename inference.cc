/*  inference.cc - General inference technique base class

 Adrian Groves & Michael Chappell, FMRIB Image Analysis Group & IBME QuBIc Group

 Copyright (C) 2007-20015 University of Oxford  */

/*  CCOPYRIGHT */

#include "inference.h"
#include "easylog.h"

#include <newmat.h>

#include <fstream>
#include <iomanip>
#include <math.h>

using namespace std;
using namespace NEWMAT;

std::vector<std::string> InferenceTechnique::GetKnown()
{
    InferenceTechniqueFactory *factory = InferenceTechniqueFactory::GetInstance();
    return factory->GetNames();
}

InferenceTechnique *InferenceTechnique::NewFromName(const string &name)
{
    InferenceTechniqueFactory *factory = InferenceTechniqueFactory::GetInstance();
    InferenceTechnique *inf = factory->Create(name);
    if (inf == NULL)
    {
        throw InvalidOptionValue("method", name, "Unrecognized inference method");
    }
    return inf;
}

void InferenceTechnique::UsageFromName(const string &name, std::ostream &stream)
{
    stream << "Usage information for method: " << name << endl << endl;

    std::auto_ptr<InferenceTechnique> method(NewFromName(name));
    stream << method->GetDescription() << endl << endl << "Options: " << endl << endl;
    vector<OptionSpec> options;
    method->GetOptions(options);
    if (options.size() > 0)
    {
        for (vector<OptionSpec>::iterator iter = options.begin(); iter != options.end(); ++iter)
        {
            stream << *iter << endl;
        }
    }
}

InferenceTechnique::InferenceTechnique()
    : m_model(NULL)
    , m_num_params(0)
    , m_halt_bad_voxel(true)
{
}

void InferenceTechnique::Initialize(FwdModel *fwd_model, FabberRunData &rundata)
{
    m_log = rundata.GetLogger();
    m_debug = rundata.GetBool("debug");
    if (m_debug)
        LOG << setprecision(17);

    m_model = fwd_model;
    vector<Parameter> params;
    m_model->GetParameters(rundata, params);

    m_num_params = params.size();
    LOG << "InferenceTechnique::Model has " << m_num_params << " parameters" << endl;

    // Read masked time points option and log if any have been specified
    m_masked_tpoints = rundata.GetIntList("mt", 1);
    if (m_masked_tpoints.size() > 0)
    {
        LOG << "InferenceTechnique::Masking " << m_masked_tpoints.size() << " time points: ";
        for (unsigned int i = 0; i < m_masked_tpoints.size(); i++)
        {
            LOG << m_masked_tpoints[i] << " ";
        }
        LOG << endl;
    }

    // Allow calculation to continue even with bad voxels
    // Note that this is a bad idea when spatialDims>0, because the bad voxel
    // will drag its neighbours around... but can never recover!  Maybe a more
    // sensible approach is to reset bad voxels to the prior on each iteration.
    // Spatial VB ignores this parameter, presumably for the above reason
    m_halt_bad_voxel = !rundata.GetBool("allow-bad-voxels");
    if (m_halt_bad_voxel)
    {
        LOG << "InferenceTechnique::Note: numerical errors in voxels will cause the program to "
               "halt.\n"
            << "InferenceTechnique::Use --allow-bad-voxels (with caution!) to keep on "
               "calculating.\n";
    }
    else
    {
        LOG << "InferenceTechnique::Using --allow-bad-voxels: numerical errors in a voxel will\n"
            << "InferenceTechnique::simply stop the calculation of that voxel.\n"
            << "InferenceTechnique::Check log for 'Going on to the next voxel' messages.\n"
            << "InferenceTechnique::Note that you should get very few (if any) exceptions like "
               "this;"
            << "InferenceTechnique::they are probably due to bugs or a numerically unstable model.";
    }
}

void InferenceTechnique::SaveResults(FabberRunData &rundata) const
{
    LOG << "InferenceTechnique::Preparing to save results..." << endl;
    int nVoxels = resultMVNs.size();

    // Save the resultMVNs NIFTI file
    if (rundata.GetBool("save-mvn"))
    {
        MVNDist::Save(resultMVNs, "finalMVN", rundata);
    }

    // Create individual files for each parameter's mean and Z-stat
    vector<Parameter> params;
    m_model->GetParameters(rundata, params);
    if (rundata.GetBool("save-mean") | rundata.GetBool("save-std") | rundata.GetBool("save-zstat") | rundata.GetBool("save-var"))
    {
        LOG << "InferenceTechnique::Writing means..." << endl;
        for (unsigned i = 1; i <= params.size(); i++)
        {
            Matrix paramMean, paramZstat, paramStd, paramVar;
            paramMean.ReSize(1, nVoxels);
            paramZstat.ReSize(1, nVoxels);
            paramStd.ReSize(1, nVoxels);
            paramVar.ReSize(1, nVoxels);

            for (int vox = 1; vox <= nVoxels; vox++)
            {
                MVNDist result = *resultMVNs[vox - 1];
                m_model->ToModel(result);
                paramMean(1, vox) = result.means(i);
                paramVar(1, vox) = result.GetCovariance()(i, i);
                double std = sqrt(paramVar(1, vox));
                paramZstat(1, vox) = paramMean(1, vox) / std;
                paramStd(1, vox) = std;
            }

            if (rundata.GetBool("save-mean"))
                rundata.SaveVoxelData("mean_" + params.at(i - 1).name, paramMean);
            if (rundata.GetBool("save-zstat"))
                rundata.SaveVoxelData("zstat_" + params.at(i - 1).name, paramZstat);
            if (rundata.GetBool("save-std"))
                rundata.SaveVoxelData("std_" + params.at(i - 1).name, paramStd);
            if (rundata.GetBool("save-var"))
                rundata.SaveVoxelData("var_" + params.at(i - 1).name, paramVar);
        }
    }

    // Produce the model fit and residual volume series
    bool saveModelFit = rundata.GetBool("save-model-fit");
    bool saveResiduals = rundata.GetBool("save-residuals");
    vector<string> outputs;
    outputs.push_back("");
    m_model->GetOutputs(outputs);
    if (saveModelFit || saveResiduals || (outputs.size() > 1))
    {
        LOG << "InferenceTechnique::Writing model time series data (fit, residuals and "
               "model-specific output)"
            << endl;

        Matrix result, residuals, datamtx, coords, suppdata;
        datamtx = rundata.GetMainVoxelData(); // it is just possible that the model needs the data
                                              // in its calculations
        coords = rundata.GetVoxelCoords();
        suppdata = rundata.GetVoxelSuppData();
        result.ReSize(datamtx.Nrows(), nVoxels);
        ColumnVector tmp;
        for (vector<string>::iterator iter = outputs.begin(); iter != outputs.end(); ++iter)
        {
            LOG << "InferenceTechnique::Evaluating model for output: " << *iter << endl;
            for (int vox = 1; vox <= nVoxels; vox++)
            {
                // do the evaluation
                try
                {
                    // pass in stuff that the model might need
                    ColumnVector y = datamtx.Column(vox);
                    ColumnVector vcoords = coords.Column(vox);
                    if (suppdata.Ncols() > 0)
                    {
                        m_model->PassData(vox, y, vcoords, suppdata.Column(vox));
                    }
                    else
                    {
                        m_model->PassData(vox, y, vcoords);
                    }

                    m_model->EvaluateFabber(
                        resultMVNs.at(vox - 1)->means.Rows(1, m_num_params), tmp, *iter);
                    if (result.Nrows() != tmp.Nrows())
                    {
                        // Only occurs on first voxel if output size is not equal to
                        // data size
                        result.ReSize(tmp.Nrows(), nVoxels);
                    }
                    result.Column(vox) = tmp;
                }
                // Ignore exceptions for the default Evaluate key - errors when evaluating the model would already have
                // occurred during inference and the relevant warnings output.
                catch (NEWMAT::Exception &e)
                {
                    if (*iter != "")
                    {
                        LOG << "InferenceTechnique::NEWMAT error generating output " << *iter << " for voxel " << vox << " : " << e.what() << endl;
                    }
                }
                catch (std::exception &e)
                {
                    if (*iter != "")
                    {
                        LOG << "InferenceTechnique::Error generating output " << *iter << " for voxel " << vox << " : " << e.what() << endl;
                    }
                }
                catch (...)
                {
                    if (*iter != "")
                    {
                        LOG << "InferenceTechnique::Unexpected error generating output " << *iter << " for voxel " << vox << " : no message available" << endl;
                    }
                }
            }
            if (*iter == "")
            {
                if (saveResiduals)
                {
                    LOG << "InferenceTechnique::Saving residuals" << endl;
                    residuals = datamtx - result;
                    rundata.SaveVoxelData("residuals", residuals);
                }
                if (saveModelFit)
                {
                    LOG << "InferenceTechnique::Saving model fit" << endl;
                    rundata.SaveVoxelData("modelfit", result);
                }
            }
            else if (rundata.GetBool("save-model-extras"))
            {
                LOG << "InferenceTechnique::Saving extra output (size=" << result.Nrows() << ")" << endl;
                rundata.SaveVoxelData(*iter, result);
            }
        }
    }

#if 0
	{
		LOG << "InferenceTechnique::Writing model variances..." << endl;
		Matrix modelStd;
		Matrix datamtx = data.GetMainVoxelData();
		modelStd.ReSize(datamtx.Nrows(), nVoxels);
		for (int vox = 1; vox <= nVoxels; vox++)
		{
			LinearizedFwdModel lin(model);
			lin.ReCentre(resultMVNs.at(vox - 1)->means.Rows(1, m_num_params));
			Matrix var = resultMVNs.at(vox - 1)->GetCovariance().SymSubMatrix(1, m_num_params);

			Matrix mvar = lin.Jacobian() * var * lin.Jacobian().t();
			ColumnVector tmp(datamtx.Nrows());
			for (int i = 1; i <= datamtx.Nrows(); i++)
			{
				tmp(i) = sqrt(mvar(i, i));
			}

			modelStd.Column(vox) = tmp;
		}
		data.SaveVoxelData("modelstd", modelStd);

	}
#endif

    LOG << "InferenceTechnique::Done writing results." << endl;
}

void InferenceTechnique::InitMVNFromFile(FabberRunData &rundata, string paramFilename = "")
{
    // Loads in a MVN to set it as inital values for inference
    // can cope with the special scenario in which extra parameters have been added to the inference
    LOG << "InferenceTechnique::Merging supplied MVN with model intialization." << endl;

    if (paramFilename == "")
    {
        MVNDist::Load(resultMVNs, "continue-from-mvn", rundata, m_log);
    }
    else
    {
        // load in parameters FIXME this will not work right now
        LOG << "InferenceTechnique::Parameters named in file" << endl;

        // Get a list of parameter names which were associated with the previous
        // run so we can merge the MVNs with the full parameter set
        string currparam;
        ifstream paramFile((paramFilename).c_str());
        if (!paramFile.good())
        {
            throw InvalidOptionValue("", paramFilename, "Could not open parameter file");
        }
        vector<string> paramNames;
        LOG << "InferenceTechnique::Parameters from previous run: " << endl;
        while (paramFile.good())
        {
            getline(paramFile, currparam);
            paramNames.push_back(currparam);
            LOG << currparam << endl;
        }
        paramNames.pop_back(); // remove final empty line assocaited with eof

        // get the parameters in the model
        vector<string> ModelparamNames;
        m_model->NameParams(ModelparamNames);
        LOG << "InferenceTechnique::Parameters named in model" << endl;
        for (int p = 0; p < m_num_params; p++)
        {
            LOG << ModelparamNames[p] << endl;
        }

        // load in the MVN
        vector<MVNDist *> MVNfile;
        MVNDist::Load(MVNfile, "continue-from-mvn", rundata, m_log);

        // Get defaults from the model. The prior is not used, the posterior is used
        // if we don't have a posterior for a parameter in the file
        MVNDist tempprior(m_num_params);
        MVNDist tempposterior(m_num_params);
        m_model->HardcodedInitialDists(tempprior, tempposterior);

        // go through the parameters in the model and either:
        // 1.) load the MVN from MVNfile if it is included, or
        // 2.) use the default value from the model

        // first work out where parameters in file MVN go in the model
        LOG << "InferenceTechnique::Matching parameters from file with model:" << endl;
        // Vector of bools for each parameter in the model. True if found in the file we're loading
        vector<bool> usefile(m_num_params, false);
        // If true, location of parameter in file (starting at 0)
        vector<int> oldloc(m_num_params, 0);
        // Vector of bools, one for each parameter in the file. True to flag matched to a file
        // parameter
        vector<bool> hasmatched(m_num_params, false);
        for (int p = 0; p < m_num_params; p++)
        {
            usefile[p] = false;
            for (unsigned q = 0; q < paramNames.size(); q++)
            {
                if (ModelparamNames[p] == paramNames[q])
                {
                    usefile[p] = true;
                    oldloc[p] = q;
                    hasmatched[q] = true;
                    LOG << ModelparamNames[p] << ": Matched with file" << endl;
                }
            }
            if (!usefile[p])
            {
                LOG << ModelparamNames[p] << ": Not matched, set from model default" << endl;
            }
        }

        // Make a note of any parameters in the file that were not matched to params in the model
        for (unsigned int q = 0; q < paramNames.size(); q++)
        {
            if (!hasmatched[q])
            {
                LOG_ERR(paramNames[q] + ": Not matched!");
            }
        }

        // for (int a=0; a<usefile.size(); a++) {
        //       cout << usefile[a] << "  " << oldloc[a] << endl;
        //     }

        ColumnVector modelmeans = tempposterior.means;
        ColumnVector newmeans = modelmeans;
        SymmetricMatrix modelcov = tempposterior.GetCovariance();
        SymmetricMatrix newcov = modelcov;

        // Number of forward model and noise params in the FILE we are loading
        int n_file_params = paramNames.size();
        int n_file_noiseparams = MVNfile[1]->means.Nrows() - n_file_params;
        int nvox = MVNfile.size();

        // For every voxel we need to add an MVNDist to resultMVNs.
        for (int v = 0; v < nvox; v++)
        {
            // Get the file data
            MVNDist fwddist = MVNfile[v]->GetSubmatrix(1, n_file_params);
            MVNDist noisedist
                = MVNfile[v]->GetSubmatrix(n_file_params + 1, n_file_params + n_file_noiseparams);

            for (unsigned int p = 0; p < ModelparamNames.size(); p++)
            {
                // deal with the means
                if (usefile[p])
                {
                    newmeans(p + 1) = fwddist.means(oldloc[p] + 1);
                }
            }
            MVNDist newfwd(m_num_params);
            newfwd.means = newmeans;

            // deal with the covariances
            SymmetricMatrix filecov = fwddist.GetCovariance();
            for (unsigned int p = 0; p < ModelparamNames.size(); p++)
            {
                for (unsigned int q = 0; q <= p; q++)
                {
                    if (usefile[p])
                    {
                        if (usefile[q])
                        {
                            newcov(p + 1, q + 1) = filecov(oldloc[p] + 1, oldloc[q] + 1);
                        }
                    }
                }
            }
            newfwd.SetCovariance(newcov);

            // Note that this assumes we have the same number of noise parameters in the
            // file as we want in the calculation, not sure at moment if this is always
            // true or not.
            MVNDist *distout = new MVNDist(newfwd, noisedist);
            resultMVNs.push_back(distout);
        }
    }
}

InferenceTechnique::~InferenceTechnique()
{
    while (!resultMVNs.empty())
    {
        delete resultMVNs.back();
        resultMVNs.pop_back();
    }
}
