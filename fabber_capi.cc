/*  fabber_capi.cc - Pure C API for fabber

 Used by language bindings, e.g. Python interface

 Martin Craig

 Copyright (C) 2016 University of Oxford  */

/*  CCOPYRIGHT */

#include "fabber_capi.h"

#include "easylog.h"
#include "fwdmodel.h"
#include "inference.h"
#include "rundata_array.h"
#include "setup.h"

#include "armawrap/newmat.h"

#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <stdio.h>
#include <string.h>

using NEWMAT::Matrix;
using namespace std;

static int fabber_err(int code, const char *msg, char *err_buf)
{
    // Error buffer is optional
    if (!err_buf)
        return code;
    if (!msg)
        msg = "NULL message";

    strncpy(err_buf, msg, FABBER_ERR_MAXC - 1);
    err_buf[FABBER_ERR_MAXC - 1] = '\0';
    return code;
}

void *fabber_new(char *err_buf)
{
    try
    {
        FabberSetup::SetupDefaults();
        FabberRunDataArray *rundata = new FabberRunDataArray(false);
        return rundata;
    }
    catch (...)
    {
        fabber_err(FABBER_ERR_FATAL, "Failed to allocate memory for run data", err_buf);
        return NULL;
    }
}

int fabber_load_models(void *fab, const char *libpath, char *err_buf)
{
    if (!fab)
        return fabber_err(FABBER_ERR_FATAL, "Rundata is NULL", err_buf);
    if (!libpath)
        return fabber_err(FABBER_ERR_FATAL, "Library path is NULL", err_buf);

    try
    {
        FwdModel::LoadFromDynamicLibrary(libpath);
        return 0;
    }
    catch (exception &e)
    {
        return fabber_err(FABBER_ERR_FATAL, e.what(), err_buf);
    }
}

int fabber_set_extent(
    void *fab, unsigned int nx, unsigned int ny, unsigned int nz, const int *mask, char *err_buf)
{
    if (!fab)
        return fabber_err(FABBER_ERR_FATAL, "Rundata is NULL", err_buf);
    if (!mask)
        return fabber_err(FABBER_ERR_FATAL, "Mask is NULL", err_buf);
    if ((nx <= 0) || (ny <= 0) || (nz <= 0))
        return fabber_err(FABBER_ERR_FATAL, "Dimensions must be >0", err_buf);

    try
    {
        FabberRunDataArray *rundata = (FabberRunDataArray *)fab;
        rundata->SetExtent(nx, ny, nz, mask);
        return 0;
    }
    catch (NEWMAT::Exception &e)
    {
        return fabber_err(FABBER_ERR_FATAL, e.what(), err_buf);
    }
    catch (exception &e)
    {
        return fabber_err(FABBER_ERR_FATAL, e.what(), err_buf);
    }
    catch (...)
    {
        return fabber_err(FABBER_ERR_FATAL, "Error setting extent", err_buf);
    }
}

int fabber_set_opt(void *fab, const char *key, const char *value, char *err_buf)
{
    if (!fab)
        return fabber_err(FABBER_ERR_FATAL, "Rundata is NULL", err_buf);

    if (!key || !value)
        return fabber_err(FABBER_ERR_FATAL, "Option key or value is NULL", err_buf);

    try
    {
        FabberRunDataArray *rundata = (FabberRunDataArray *)fab;
        rundata->Set(key, value);
        return 0;
    }
    catch (exception &e)
    {
        return fabber_err(FABBER_ERR_FATAL, e.what(), err_buf);
    }
}

int fabber_set_data(
    void *fab, const char *name, unsigned int data_size, const float *data, char *err_buf)
{
    if (!fab)
        return fabber_err(FABBER_ERR_FATAL, "Rundata is NULL", err_buf);
    if (!data)
        return fabber_err(FABBER_ERR_FATAL, "Data buffer is NULL", err_buf);
    if (!name)
        return fabber_err(FABBER_ERR_FATAL, "Data name is NULL", err_buf);
    if (data_size <= 0)
        return fabber_err(FABBER_ERR_FATAL, "Data size must be >0", err_buf);

    FabberRunDataArray *rundata = (FabberRunDataArray *)fab;
    EasyLog log;
    rundata->SetLogger(&log); // Ignored but avoids it going to stdout
    try
    {
        rundata->SetVoxelDataArray(name, data_size, data);
    }
    catch (exception &e)
    {
        return fabber_err(FABBER_ERR_FATAL, e.what(), err_buf);
    }
    catch (...)
    {
        return fabber_err(FABBER_ERR_FATAL, "Error setting data", err_buf);
    }
    return 0;
}

int fabber_get_data_size(void *fab, const char *name, char *err_buf)
{
    if (!fab)
        return fabber_err(FABBER_ERR_FATAL, "Rundata is NULL", err_buf);
    if (!name)
        return fabber_err(FABBER_ERR_FATAL, "Data name is NULL", err_buf);

    FabberRunDataArray *rundata = (FabberRunDataArray *)fab;
    EasyLog log;
    rundata->SetLogger(&log); // Ignored but avoids it going to stdout
    try
    {
        return rundata->GetVoxelDataSize(name);
    }
    catch (DataNotFound &e)
    {
        return fabber_err(-1, e.what(), err_buf);
    }
    catch (...)
    {
        return fabber_err(FABBER_ERR_FATAL, "Error getting data", err_buf);
    }
}

int fabber_get_data(void *fab, const char *name, float *data_buf, char *err_buf)
{
    if (!fab)
        return fabber_err(FABBER_ERR_FATAL, "Rundata is NULL", err_buf);
    if (!name)
        return fabber_err(FABBER_ERR_FATAL, "Data name is NULL", err_buf);
    if (!data_buf)
        return fabber_err(FABBER_ERR_FATAL, "Data name is NULL", err_buf);

    FabberRunDataArray *rundata = (FabberRunDataArray *)fab;
    EasyLog log;
    rundata->SetLogger(&log); // Ignored but avoids it going to stdout
    try
    {
        rundata->GetVoxelDataArray(name, data_buf);
        return 0;
    }
    catch (DataNotFound &e)
    {
        return fabber_err(-1, e.what(), err_buf);
    }
    catch (...)
    {
        return fabber_err(FABBER_ERR_FATAL, "Error getting data", err_buf);
    }
}

int fabber_dorun(void *fab, unsigned int log_bufsize, char *log_buf, char *err_buf,
    void (*progress_cb)(int, int))
{
    EasyLog log;

    if (!fab)
        return fabber_err(FABBER_ERR_FATAL, "Rundata is NULL", err_buf);
    if (!log_buf)
        return fabber_err(FABBER_ERR_FATAL, "Log buffer is NULL", err_buf);
    if (!err_buf)
        return fabber_err(FABBER_ERR_FATAL, "Error buffer is NULL", err_buf);
    if (log_bufsize > 0 && !log_buf)
        return fabber_err(FABBER_ERR_FATAL, "Log buffer is NULL", err_buf);

    int ret = 0;
    FabberRunDataArray *rundata = (FabberRunDataArray *)fab;
    rundata->SetLogger(&log);
    stringstream logstr;
    try
    {
        log.StartLog(logstr);
        if (progress_cb)
        {
            CallbackProgressCheck prog(progress_cb);
            rundata->Run(&prog);
        }
        else
        {
            rundata->Run();
        }
        log.ReissueWarnings();
    }
    catch (const FabberError &e)
    {
        log.ReissueWarnings();
        log.LogStream() << e.what() << endl;
        ret = fabber_err(FABBER_ERR_FATAL, e.what(), err_buf);
    }
    catch (NEWMAT::Exception &e)
    {
        log.ReissueWarnings();
        log.LogStream() << "NEWMAT exception caught in fabber:\n  " << e.what() << endl;
        ret = fabber_err(FABBER_ERR_NEWMAT, e.what(), err_buf);
    }
    catch (const exception &e)
    {
        log.ReissueWarnings();
        log.LogStream() << "STL exception caught in fabber:\n  " << e.what() << endl;
        ret = fabber_err(FABBER_ERR_FATAL, e.what(), err_buf);
    }
    catch (...)
    {
        log.ReissueWarnings();
        log.LogStream() << "Some other exception caught in fabber!" << endl;
        ret = fabber_err(FABBER_ERR_FATAL, "Unrecognized exception", err_buf);
    }

    log.StopLog();

    strncpy(log_buf, logstr.str().c_str(), log_bufsize - 1);
    log_buf[log_bufsize - 1] = '\0';

    return ret;
}

void fabber_destroy(void *fab)
{
    if (fab)
    {
        // Get rid of registered models etc
        FabberSetup::Destroy();
        FabberRunDataArray *rundata = (FabberRunDataArray *)fab;
        delete rundata;
    }
}

int fabber_get_options(void *fab, const char *key, const char *value, unsigned int out_bufsize,
    char *out_buf, char *err_buf)
{
    if (!fab)
        return fabber_err(FABBER_ERR_FATAL, "Rundata is NULL", err_buf);
    if (!out_buf)
        return fabber_err(FABBER_ERR_FATAL, "Output buffer is NULL", err_buf);
    if (key && !value)
        return fabber_err(FABBER_ERR_FATAL, "Key specified but no value", err_buf);

    vector<OptionSpec> options;
    string desc;
    try
    {
        if (!key || (strlen(key) == 0))
        {
            FabberRunDataArray::GetOptions(options);
        }
        else if (strcmp(key, "model") == 0)
        {
            std::auto_ptr<FwdModel> model(FwdModel::NewFromName(value));
            desc = model->GetDescription();
            model->GetOptions(options);
        }
        else if (strcmp(key, "method") == 0)
        {
            std::auto_ptr<InferenceTechnique> method(InferenceTechnique::NewFromName(value));
            desc = method->GetDescription();
            method->GetOptions(options);
        }

        // Remove newlines from the description so we can guarantee that the first line is the
        // description
        desc.erase(std::remove(desc.begin(), desc.end(), '\n'), desc.end());
        stringstream out;
        out << desc << endl;

        vector<OptionSpec>::iterator iter;
        for (iter = options.begin(); iter != options.end(); ++iter)
        {
            out << iter->name << "\t" << iter->description << "\t" << iter->type << "\t"
                << iter->optional << "\t" << iter->def << endl;
        }
        string outstr = out.str();

        if (outstr.size() >= out_bufsize)
        {
            return fabber_err(-1, "Buffer too small", err_buf);
        }
        strncpy(out_buf, outstr.c_str(), outstr.size());
        out_buf[outstr.size()] = '\0';
        return 0;
    }
    catch (exception &e)
    {
        return fabber_err(FABBER_ERR_FATAL, e.what(), err_buf);
    }
    catch (...)
    {
        return fabber_err(FABBER_ERR_FATAL, "Error in get_options", err_buf);
    }
}

int fabber_get_models(void *fab, unsigned int out_bufsize, char *out_buf, char *err_buf)
{
    if (!fab)
        return fabber_err(FABBER_ERR_FATAL, "Rundata is NULL", err_buf);
    if (!out_buf)
        return fabber_err(FABBER_ERR_FATAL, "Output buffer is NULL", err_buf);

    try
    {
        vector<string> known = FwdModel::GetKnown();
        stringstream out;
        vector<string>::iterator iter;
        for (iter = known.begin(); iter != known.end(); ++iter)
        {
            out << *iter << endl;
        }
        string outstr = out.str();
        if (outstr.size() >= out_bufsize)
        {
            return fabber_err(-1, "Buffer too small", err_buf);
        }
        strncpy(out_buf, outstr.c_str(), outstr.size());
        out_buf[outstr.size()] = '\0';
        return 0;
    }
    catch (exception &e)
    {
        return fabber_err(FABBER_ERR_FATAL, e.what(), err_buf);
    }
    catch (...)
    {
        return fabber_err(FABBER_ERR_FATAL, "Error in get_models", err_buf);
    }
}

int fabber_get_methods(void *fab, unsigned int out_bufsize, char *out_buf, char *err_buf)
{
    if (!fab)
        return fabber_err(FABBER_ERR_FATAL, "Rundata is NULL", err_buf);
    if (!out_buf)
        return fabber_err(FABBER_ERR_FATAL, "Output buffer is NULL", err_buf);

    try
    {
        vector<string> known = InferenceTechnique::GetKnown();
        stringstream out;
        vector<string>::iterator iter;
        for (iter = known.begin(); iter != known.end(); ++iter)
        {
            out << *iter << endl;
        }
        string outstr = out.str();
        if (outstr.size() >= out_bufsize)
        {
            return fabber_err(-1, "Buffer too small", err_buf);
        }
        strncpy(out_buf, outstr.c_str(), outstr.size());
        out_buf[outstr.size()] = '\0';
        return 0;
    }
    catch (exception &e)
    {
        return fabber_err(FABBER_ERR_FATAL, e.what(), err_buf);
    }
    catch (...)
    {
        return fabber_err(FABBER_ERR_FATAL, "Error in get_methods", err_buf);
    }
}

int fabber_get_model_params(void *fab, unsigned int out_bufsize, char *out_buf, char *err_buf)
{
    if (!fab)
        return fabber_err(FABBER_ERR_FATAL, "Rundata is NULL", err_buf);
    if (!out_buf)
        return fabber_err(FABBER_ERR_FATAL, "Output buffer is NULL", err_buf);

    try
    {
        FabberRunDataArray *rundata = (FabberRunDataArray *)fab;
        std::auto_ptr<FwdModel> model(FwdModel::NewFromName(rundata->GetString("model")));
        EasyLog log;
        model->SetLogger(&log); // We ignore the log but this stops it going to cerr
        model->Initialize(*rundata);
        vector<Parameter> params;
        model->GetParameters(*rundata, params);
        stringstream out;
        vector<Parameter>::iterator iter;
        for (iter = params.begin(); iter != params.end(); ++iter)
        {
            out << iter->name << endl;
        }
        string outstr = out.str();
        if (outstr.size() >= out_bufsize)
        {
            return fabber_err(-1, "Buffer too small", err_buf);
        }
        strncpy(out_buf, outstr.c_str(), outstr.size());
        out_buf[outstr.size()] = '\0';
        return 0;
    }
    catch (exception &e)
    {
        return fabber_err(FABBER_ERR_FATAL, e.what(), err_buf);
    }
    catch (...)
    {
        return fabber_err(FABBER_ERR_FATAL, "Error in get_model_params", err_buf);
    }
}

int fabber_get_model_param_descs(void *fab, unsigned int out_bufsize, char *out_buf, char *err_buf)
{
    if (!fab)
        return fabber_err(FABBER_ERR_FATAL, "Rundata is NULL", err_buf);
    if (!out_buf)
        return fabber_err(FABBER_ERR_FATAL, "Output buffer is NULL", err_buf);

    try
    {
        FabberRunDataArray *rundata = (FabberRunDataArray *)fab;
        std::auto_ptr<FwdModel> model(FwdModel::NewFromName(rundata->GetString("model")));
        EasyLog log;
        model->SetLogger(&log); // We ignore the log but this stops it going to cerr
        model->Initialize(*rundata);
        vector<Parameter> params;
        model->GetParameters(*rundata, params);
        stringstream out;
        vector<Parameter>::iterator iter;
        for (iter = params.begin(); iter != params.end(); ++iter)
        {
            out << iter->name << " " << iter->desc;
            if (iter->units != "")
            {
                out << " (units: " << iter->units << ")";
            }
            out << endl;
        }
        string outstr = out.str();
        if (outstr.size() >= out_bufsize)
        {
            return fabber_err(-1, "Buffer too small", err_buf);
        }
        strncpy(out_buf, outstr.c_str(), outstr.size());
        out_buf[outstr.size()] = '\0';
        return 0;
    }
    catch (exception &e)
    {
        return fabber_err(FABBER_ERR_FATAL, e.what(), err_buf);
    }
    catch (...)
    {
        return fabber_err(FABBER_ERR_FATAL, "Error in get_model_params", err_buf);
    }
}

int fabber_get_model_outputs(void *fab, unsigned int out_bufsize, char *out_buf, char *err_buf)
{
    if (!fab)
        return fabber_err(FABBER_ERR_FATAL, "Rundata is NULL", err_buf);
    if (!out_buf)
        return fabber_err(FABBER_ERR_FATAL, "Output buffer is NULL", err_buf);

    try
    {
        FabberRunDataArray *rundata = (FabberRunDataArray *)fab;
        std::auto_ptr<FwdModel> model(FwdModel::NewFromName(rundata->GetString("model")));
        EasyLog log;
        model->SetLogger(&log); // We ignore the log but this stops it going to cerr
        model->Initialize(*rundata);
        vector<string> outputs;
        model->GetOutputs(outputs);
        stringstream out;
        vector<string>::iterator iter;
        for (iter = outputs.begin(); iter != outputs.end(); ++iter)
        {
            if (*iter != "")
                out << *iter << endl;
        }
        string outstr = out.str();
        if (outstr.size() >= out_bufsize)
        {
            return fabber_err(-1, "Buffer too small", err_buf);
        }
        strncpy(out_buf, outstr.c_str(), outstr.size());
        out_buf[outstr.size()] = '\0';
        return 0;
    }
    catch (exception &e)
    {
        return fabber_err(FABBER_ERR_FATAL, e.what(), err_buf);
    }
    catch (...)
    {
        return fabber_err(FABBER_ERR_FATAL, "Error in fabber_get_model_outputs", err_buf);
    }
}

int fabber_model_evaluate(void *fab, unsigned int n_params, float *params, unsigned int n_ts,
    float *indata, float *output, char *err_buf)
{
    return fabber_model_evaluate_output(fab, n_params, params, n_ts, indata, "", output, err_buf);
}

int fabber_model_evaluate_output(void *fab, unsigned int n_params, float *params, unsigned int n_ts,
    float *indata, const char *output_name, float *output, char *err_buf)
{
    if (!fab)
        return fabber_err(FABBER_ERR_FATAL, "Rundata is NULL", err_buf);
    if (!params)
        return fabber_err(FABBER_ERR_FATAL, "Params array is NULL", err_buf);
    if (!output_name)
        return fabber_err(FABBER_ERR_FATAL, "Output name is NULL", err_buf);
    if (!output)
        return fabber_err(FABBER_ERR_FATAL, "Output buffer is NULL", err_buf);

    EasyLog log;
    int ret = FABBER_ERR_FATAL;
    stringstream logstr;
    try
    {
        FabberRunDataArray *rundata = (FabberRunDataArray *)fab;
        std::auto_ptr<FwdModel> model(FwdModel::NewFromName(rundata->GetString("model")));

        log.StartLog(logstr);
        model->SetLogger(&log);
        model->Initialize(*rundata);
        NEWMAT::ColumnVector p_vec(n_params);
        NEWMAT::ColumnVector o_vec(n_ts);
        NEWMAT::ColumnVector data_vec(n_ts);
        NEWMAT::ColumnVector coords(3);
        for (unsigned int i = 0; i < n_params; i++)
        {
            p_vec(i + 1) = params[i];
        }
        for (unsigned int i=0; i < n_ts; i++)
        {
            if (indata)
                data_vec(i + 1) = indata[i];
            else
                data_vec(i + 1) = 0;
        }
        coords(1) = 1;
        coords(2) = 1;
        coords(3) = 1;
        model->PassData(1, data_vec, coords);
        model->EvaluateModel(p_vec, o_vec, output_name);
        for (unsigned int i = 0; i < n_ts; i++)
        {
            // Model may not return the same number of timepoints as passed in!
            if ((int)i < o_vec.Nrows())
                output[i] = o_vec(i + 1);
            else
                output[i] = 0;
        }

        log.ReissueWarnings();
        ret = 0;
    }
    catch (NEWMAT::Exception &e)
    {
        ret = fabber_err(FABBER_ERR_FATAL, e.what(), err_buf);
    }
    catch (exception &e)
    {
        ret = fabber_err(FABBER_ERR_FATAL, e.what(), err_buf);
    }
    catch (...)
    {
        ret = fabber_err(FABBER_ERR_FATAL, "Error evaluating model", err_buf);
    }

    log.StopLog();
    return ret;
}
