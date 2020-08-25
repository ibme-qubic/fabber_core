/*  dataset.cc - Data-loading class for FABBER

 Adrian Groves, FMRIB Image Analysis Group

 Copyright (C) 2007-2008 University of Oxford  */

/*  CCOPYRIGHT */

#include "rundata.h"

#include "easylog.h"
#include "fwdmodel.h"
#include "inference.h"
#include "setup.h"
#include "version.h"

#include "armawrap/newmat.h"

#include <errno.h>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef _WIN32
#include "direct.h"
#endif

#include <sys/stat.h>
#include <time.h>

using namespace std;
using namespace NEWMAT;

static bool is_dir(string path)
{
#ifdef _WIN32
    struct _stat s;
    if (_stat(path.c_str(), &s) == 0)
#else
    struct stat s;
    if (stat(path.c_str(), &s) == 0)
#endif
    {
        return (s.st_mode & S_IFDIR);
    }
    else
    {
        // Does not exist, so not a directory...
        return false;
    }
}

std::ostream &operator<<(std::ostream &out, const OptionType value)
{
    const char *s = 0;
    switch (value)
    {
    case OPT_BOOL:
        s = "BOOL";
        break;
    case OPT_STR:
        s = "STR";
        break;
    case OPT_INT:
        s = "INT";
        break;
    case OPT_FLOAT:
        s = "FLOAT";
        break;
    case OPT_FILE:
        s = "FILE";
        break;
    case OPT_IMAGE:
        s = "IMAGE";
        break;
    case OPT_TIMESERIES:
        s = "TIMESERIES";
        break;
    case OPT_MVN:
        s = "MVN";
        break;
    case OPT_MATRIX:
        s = "MATRIX";
        break;
    default:
        s = "UNKNOWN";
    }

    return out << s;
}

std::ostream &operator<<(std::ostream &out, const OptionSpec &value)
{
    return out << "--" << value.name << " [" << value.type << ","
               << (value.optional ? "NOT REQUIRED" : "REQUIRED") << ","
               << ((value.def == "") ? "NO DEFAULT" : "DEFAULT=" + value.def) << "]" << endl
               << "        " << value.description << endl;
}

void PercentProgressCheck::Progress(int voxel, int nVoxels)
{
    if (nVoxels == 0)
    {
        cout << "100%" << endl;
        return;
    }

    int percent = (100 * voxel) / nVoxels;
    if (percent > m_last)
    {
        cout << "\b\b\b";
        m_last = percent;
        if (m_last == 0)
            cout << " ";
        cout << m_last << "%" << flush;
        if (m_last == 100)
            cout << endl;
    }
}

void SimpleProgressCheck::Progress(int voxel, int nVoxels)
{
    if (nVoxels == 0)
    {
        cout << "100" << endl;
        return;
    }

    int percent = (100 * voxel) / nVoxels;
    if (percent > m_last)
    {
        m_last = percent;
        cout << m_last << endl << flush;
    }
}

static OptionSpec OPTIONS[] = {
    { "help", OPT_BOOL, "Print this usage method. If given with --method or --model, display "
                        "relevant method/model usage information",
        OPT_NONREQ, "" },
    { "listmethods", OPT_BOOL, "List all known inference methods", OPT_NONREQ, "" },
    { "listmodels", OPT_BOOL, "List all known forward models", OPT_NONREQ, "" },
    { "listparams", OPT_BOOL, "List model parameters (requires model configuration options to be given)", OPT_NONREQ, "" },
    { "descparams", OPT_BOOL, "Descript model parameters (name, description, units) - requires model configuration options to be given. Note that not all models provide parameter descriptions", OPT_NONREQ, "" },
    { "listoutputs", OPT_BOOL, "List additional model outputs (requires model configuration options to be given)", OPT_NONREQ, "" },
    { "evaluate", OPT_STR, "Evaluate model. Set to name of output required or blank for default output. Requires model configuration options, --evaluate-params and --evaluate-nt", OPT_NONREQ, "" },
    { "evaluate-params", OPT_MATRIX, "List of parameter values for evaluation", OPT_NONREQ, "" },
    { "evaluate-nt", OPT_INT, "Number of time points for evaluation - must be consistent with model options where appropriate", OPT_NONREQ, "" },
    { "simple-output", OPT_BOOL, "Instead of usual standard output, simply output series of lines each giving progress as percentage", OPT_NONREQ, "" },
    { "output", OPT_STR, "Directory for output files (including logfile)", OPT_REQ, "" },
    { "overwrite", OPT_BOOL, "If set will overwrite existing output. If not set, new output "
                             "directories will be created by appending '+' to the directory name ",
        OPT_NONREQ, "" },
    { "link-to-latest", OPT_BOOL,
        "Try to create a link to the most recent output directory with the prefix _latest",
        OPT_NONREQ, "" },
    { "method", OPT_STR, "Use this inference method", OPT_REQ, "" },
    { "model", OPT_STR, "Use this forward model", OPT_REQ, "" },
    { "loadmodels", OPT_FILE,
        "Load models dynamically from the specified filename, which should be a DLL/shared library",
        OPT_NONREQ, "" },
    { "data", OPT_TIMESERIES, "Specify a single input data file", OPT_REQ, "" },
    { "data<n>", OPT_TIMESERIES, "Specify multiple data files for n=1, 2, 3...", OPT_NONREQ, "" },
    { "data-order", OPT_STR, "If multiple data files are specified, how they will be handled: "
                             "concatenate = one after the other,  interleave = first record from "
                             "each file, then  second, etc.",
        OPT_NONREQ, "interleave" },
    { "mask", OPT_IMAGE, "Mask file. Inference will only be performed where mask value > 0",
        OPT_NONREQ, "" },
    { "mt<n>", OPT_INT, "List of masked time points, indexed from 1. These will be ignored in the "
                        "parameter updates",
        OPT_NONREQ, "" },
    { "suppdata", OPT_TIMESERIES, "'Supplemental' timeseries data, required for some models",
        OPT_NONREQ, "" },
    { "dump-param-names", OPT_BOOL,
        "Write the file paramnames.txt containing the names of the model parameters", OPT_NONREQ,
        "" },
    { "save-model-fit", OPT_BOOL, "Output the model prediction as a 4d volume", OPT_NONREQ, "" },
    { "save-residuals", OPT_BOOL,
        "Output the residuals (difference between the data and the model prediction)", OPT_NONREQ,
        "" },
    { "save-model-extras", OPT_BOOL, "Output any additional model-specific timeseries data",
        OPT_NONREQ, "" },
    { "save-mvn", OPT_BOOL, "Output the final MVN distributions.", OPT_NONREQ, "" },
    { "save-mean", OPT_BOOL, "Output the parameter means.", OPT_NONREQ, "" },
    { "save-std", OPT_BOOL, "Output the parameter standard deviations.", OPT_NONREQ, "" },
    { "save-var", OPT_BOOL, "Output the parameter variances.", OPT_NONREQ, "" },
    { "save-zstat", OPT_BOOL, "Output the parameter Zstats.", OPT_NONREQ, "" },
    { "save-noise-mean", OPT_BOOL, "Output the noise means. The noise distribution inferred is the precision of a Gaussian noise source", OPT_NONREQ, "" },
    { "save-noise-std", OPT_BOOL, "Output the noise standard deviations. ", OPT_NONREQ, "" },
    { "save-free-energy", OPT_BOOL, "Output the free energy, if calculated. ", OPT_NONREQ, "" },
    { "optfile", OPT_BOOL, "File containing additional options, one per line, in the same form as "
                           "specified on the command line",
        OPT_NONREQ, "" },
    { "debug", OPT_BOOL,
        "Output large amounts of debug information. ONLY USE WITH VERY SMALL NUMBERS OF VOXELS",
        OPT_NONREQ, "" },
    { "" },
};

void FabberRunData::GetOptions(std::vector<OptionSpec> &opts)
{
    for (int i = 0; OPTIONS[i].name != ""; i++)
    {
        opts.push_back(OPTIONS[i]);
    }
}

FabberRunData::FabberRunData(bool compat_options)
    : m_progress(0)
{
    init(compat_options);
}

void FabberRunData::init(bool compat_options)
{
    FabberSetup::SetupDefaults();

    if (compat_options && !GetBool("no-compat-output"))
    {
        // For backwards compatibility with previous version of Fabber, save these data items by
        // default
        SetBool("save-mean");
        SetBool("save-std");
        SetBool("save-zstat");
        SetBool("save-noise-mean");
        SetBool("save-noise-std");
        SetBool("save-free-energy");
        SetBool("save-mvn");
    }
}

FabberRunData::~FabberRunData()
{
}

void FabberRunData::LogParams()
{
    map<string, string>::iterator iter;
    for (iter = m_params.begin(); iter != m_params.end(); ++iter)
    {
        LOG << "FabberRunData::Parameter " << iter->first << "=" << iter->second << endl;
    }
}

void FabberRunData::Run(ProgressCheck *progress)
{
    if (!m_log)
    {
        m_log = &m_default_log;
    }

    m_progress = progress;
    LOG << "FabberRunData::FABBER release: " << fabber_version() << endl;
    LOG << "FabberRunData::Last commit: " << fabber_source_date() << endl;

    time_t startTime;
    time(&startTime);
    LOG << "FabberRunData::Start time: " << ctime(&startTime);

    LogParams();

    // Set the forward model
    std::auto_ptr<FwdModel> fwd_model(FwdModel::NewFromName(GetString("model")));

    // For backwards compatibility - model may not have called superclass initialize
    fwd_model->SetLogger(m_log);
    fwd_model->Initialize(*this);

    std::vector<Parameter> params;
    fwd_model->GetParameters(*this, params);
    assert(params.size() > 0);
    LOG << "FabberRunData::Forward Model version " << fwd_model->ModelVersion() << endl;

    // Write the paramnames.txt file if required
    if (GetBool("dump-param-names"))
    {
        ofstream paramFile((GetStringDefault("output", ".") + "/paramnames.txt").c_str());
        for (unsigned i = 0; i < params.size(); i++)
        {
            paramFile << params[i].name << endl;
        }
        paramFile.close();
    }

    // Set the inference technique (and pass in the model)
    std::auto_ptr<InferenceTechnique> infer(InferenceTechnique::NewFromName(GetString("method")));
    infer->Initialize(fwd_model.get(), *this);

    // Calculations
    int nvoxels = GetVoxelCoords().Ncols();
    LOG << "FabberRunData::Num voxels " << nvoxels << endl;
    Progress(0, nvoxels);
    infer->DoCalculations(*this);
    Progress(nvoxels, nvoxels);
    LOG << "FabberRunData::Saving results " << endl;
    infer->SaveResults(*this);

    LOG << "FabberRunData::All done." << endl;

    // Options should all have been used by now, so complain if there's anything left.
    CheckAllOptionsUsed();

    time_t endTime;
    time(&endTime);
    LOG << "FabberRunData::Start time: " << ctime(&startTime); // Bizarrely, ctime() ends with a \n.
    LOG << "FabberRunData::End time: " << ctime(&endTime);
    LOG << "FabberRunData::Duration: " << endTime - startTime << " seconds." << endl;
}

static string trim(string const &str)
{
    if (str.empty())
        return str;

    size_t firstScan = str.find_first_not_of(' ');
    size_t first = firstScan == string::npos ? str.length() : firstScan;
    size_t last = str.find_last_not_of(' ');
    return str.substr(first, last - first + 1);
}

void FabberRunData::ParseParamFile(const string &filename)
{
    ifstream is(filename.c_str());
    if (!is.good())
    {
        throw FabberRunDataError("Couldn't read input options file:" + filename);
    }
    while (is.good())
    {
        string input;
        std::getline(is, input);
        input = trim(input);
        if (input.size() > 0)
        {
            if (input[0] == '#')
                continue;
            else
            {
                AddKeyEqualsValue(input, true);
            }
        }
    }
}

void FabberRunData::ParseOldStyleParamFile(const string &filename)
{
    ifstream is(filename.c_str());
    if (!is.good())
    {
        throw FabberRunDataError("Couldn't read input file: -@ " + filename);
    }
    char c;
    string param;
    while (is.good())
    {
        is.get(c);
        if (!isspace(c))
            param += c;
        else if (param == "")
        {
        }                                     // repeated whitespace, so do nothing
        else if (string(param, 0, 2) == "--") // we have an option
        {
            AddKeyEqualsValue(string(param, 2, string::npos));
            param = "";
        }
        else if (string(param, 0, 1) == "#") // comment
        {
            // discard this word and the rest of the line.
            param = "";
            while (is.good() && c != '\n')
                is.get(c);
        }
        else if (string(param, 0, 2) == "-@")
        {
            throw FabberRunDataError("Can only use -@ on the command line");
        }
        else
        {
            throw FabberRunDataError(
                "Invalid data '" + param + "' found in file '" + filename + "'");
        }
    }
}

void FabberRunData::Parse(int argc, char **argv)
{
    m_params[""] = argv[0];
    for (int a = 1; a < argc; a++)
    {
        if (string(argv[a]) == "-f")
        {
            if (++a < argc)
            {
                ParseParamFile(argv[a]);
            }
            else
            {
                throw InvalidOptionValue("-f", "", "No filename specified");
            }
        }
        else if (string(argv[a], 0, 2) == "--")
        {
            string key = argv[a] + 2; // skip the "--"
            AddKeyEqualsValue(key);
        }
        else if (string(argv[a]) == "-@")
        {
            // DEPRECATED - not POSIX compliant but kept for compatibility
            ParseOldStyleParamFile(argv[++a]);
        }
        else
        {
            throw FabberRunDataError("Option '" + string(argv[a]) + "' doesn't begin with --");
        }
    }
    if (HaveKey("optfile"))
    {
        // Replacement for non-POSIX compliant -@ option
        ParseOldStyleParamFile(GetString("optfile"));
    }
}

void FabberRunData::AddKeyEqualsValue(const string &exp, bool trim_comments)
{
    string::size_type eqPos = exp.find("=");
    string key = trim(string(exp, 0, eqPos));
    if (eqPos != (exp.npos))
    {
        string::size_type end = exp.npos;
        if (trim_comments)
            end = exp.find("#");
        string value = trim(exp.substr(eqPos + 1, end - (eqPos + 1)));
        if (m_params.count(key) > 0)
            throw InvalidOptionValue(key, value, "Already has a value: " + m_params[key]);

        if (key == "loadmodels")
        {
            FwdModel::LoadFromDynamicLibrary(value, m_log);
        }
        else
        {
            m_params[key] = value;
        }
    }
    else
    {
        m_params[exp] = "";
    }
}

void FabberRunData::Set(const string &key, const string &value)
{
    m_params[key] = value;
}
void FabberRunData::Set(const string &key, double value)
{
    m_params[key] = stringify(value);
}
void FabberRunData::SetBool(const string &key, bool value)
{
    if (value)
        m_params[key] = "";
    else
        m_params.erase(key);
}

void FabberRunData::Unset(const std::string &key)
{
    m_params.erase(key);
}
bool FabberRunData::HaveKey(const string &key)
{
    return m_params.count(key) > 0;
}
string FabberRunData::GetString(const string &key)
{
    return Read(key, key);
}
string FabberRunData::GetStringDefault(const string &key, const string &def) const
{
    // Flag parameter as having been read
    m_used_params.insert(key);
    if (m_params.count(key) == 0)
    {
        return def;
    }
    else
    {
        return m_params.find(key)->second;
    }
}

std::vector<std::string> FabberRunData::GetStringList(const std::string &prefix)
{
    std::vector<std::string> ret;
    if (HaveKey(prefix))
    {
        ret.push_back(GetString(prefix));
    }
    else
    {
        int n = 1;
        while (HaveKey(prefix + stringify(n)))
        {
            ret.push_back(GetString(prefix + stringify(n)));
            n++;
        }
    }
    return ret;
}

bool FabberRunData::GetBool(const string &key)
{
    if (m_params.count(key) == 0)
        return false;

    // Flag parameter as having been read
    m_used_params.insert(key);
    if (m_params[key] == "")
    {
        return true;
    }

    throw InvalidOptionValue(key, m_params[key], "Value should not be given for boolean option");
}

int FabberRunData::GetInt(const string &key, int min, int max)
{
    string val = GetString(key);
    try
    {
        int i = convertTo<int>(val, key);
        if (i < min)
            throw InvalidOptionValue(key, val, "Minimum " + stringify(min));
        if (i > max)
            throw InvalidOptionValue(key, val, "Maximum " + stringify(max));
        return i;
    }
    catch (invalid_argument &)
    {
        throw InvalidOptionValue(key, val, "Must be an integer");
    }
}

int FabberRunData::GetIntDefault(const string &key, int def, int min, int max)
{
    if (m_params.count(key) == 0)
        return def;
    else
        return GetInt(key, min, max);
}

std::vector<int> FabberRunData::GetIntList(const std::string &prefix, int min, int max)
{
    std::vector<int> ret;
    if (HaveKey(prefix))
    {
        ret.push_back(GetInt(prefix, min, max));
    }
    else
    {
        int n = 1;
        while (HaveKey(prefix + stringify(n)))
        {
            ret.push_back(GetInt(prefix + stringify(n), min, max));
            n++;
        }
    }
    return ret;
}

double FabberRunData::GetDouble(const string &key, double min, double max)
{
    string val = GetString(key);
    try
    {
        double d = convertTo<double>(val, key);
        if (d < min)
            throw InvalidOptionValue(key, val, "Minimum " + stringify(min));
        if (d > max)
            throw InvalidOptionValue(key, val, "Maximum " + stringify(max));
        return d;
    }
    catch (invalid_argument &)
    {
        throw InvalidOptionValue(key, val, "Must be an number");
    }
}

double FabberRunData::GetDoubleDefault(const string &key, double def, double min, double max)
{
    if (m_params.count(key) == 0)
        return def;
    else
        return GetDouble(key, min, max);
}

std::vector<double> FabberRunData::GetDoubleList(const std::string &prefix, double min, double max)
{
    std::vector<double> ret;
    if (HaveKey(prefix))
    {
        ret.push_back(GetDouble(prefix, min, max));
    }
    else
    {
        int n = 1;
        while (HaveKey(prefix + stringify(n)))
        {
            ret.push_back(GetDouble(prefix + stringify(n), min, max));
            n++;
        }
    }
    return ret;
}

string FabberRunData::Read(const string &key, const string &msg)
{
    if (m_params.count(key) == 0)
        throw MandatoryOptionMissing(msg);

    if (m_params[key] == "")
        throw InvalidOptionValue(key, "<no value>", "Value must be given");

    // Flag option as having been used and return it
    m_used_params.insert(key);
    return m_params[key];
}

std::string FabberRunData::Read(const std::string &key)
{
    return GetString(key);
}
std::string FabberRunData::ReadWithDefault(const std::string &key, const std::string &def)
{
    return GetStringDefault(key, def);
}

bool FabberRunData::ReadBool(const std::string &key)
{
    return GetBool(key);
}

void FabberRunData::CheckAllOptionsUsed() const
{
    for (map<string, string>::const_iterator iter = m_params.begin(); iter != m_params.end();
         ++iter)
    {
        if ((iter->first != "") && m_used_params.find(iter->first) == m_used_params.end())
        {
            WARN_ONCE("Unused option specified: " + iter->first);
        }
    }
}

string FabberRunData::GetOutputDir()
{
    // Avoid warning about unused parameter, this is set by default but not always relevant
    bool ltl = GetBool("link-to-latest");

    if (m_outdir != "")
        return m_outdir;

    string basename = GetStringDefault("output", "");
    if (basename == "")
    {
        m_outdir = ".";
        return m_outdir;
    }
    bool overwrite = GetBool("overwrite");

    // From Wooly's utils/log.cc
    int count = 0;
    m_outdir = basename;
    while (true)
    {
        if (count >= 50) // I'm using a lot for some things
        {
            throw FabberInternalError(
                ("Cannot create output directory (bad path, or too many + signs?): " + m_outdir)
                    .c_str());
        }

        // Clear errno so it can be inspected later; result is only meaningful if mkdir fails.
        errno = 0;
        int ret = 0;
#ifdef _WIN32
        ret = _mkdir(m_outdir.c_str());
#else
        ret = mkdir(m_outdir.c_str(), 0777);
#endif

        if (ret == 0)
        {
            // Success, directory created
            break;
        }
        else if (overwrite)
        {
            if ((errno == EEXIST) && is_dir(m_outdir))
            {
                // If directory already exists -- that's fine.
                break;
            }
            else
            {
                // Other error -- might be a problem!
                throw FabberInternalError(
                    ("Unexpected problem creating output directory in overwrite mode: " + m_outdir)
                        .c_str());
            }
        }
        m_outdir += "+";
        count++;
    }

#ifdef _WIN32
#else
    // Might be useful for jobs running on the queue:
    int ret = system(("uname -a > " + m_outdir + "/uname.txt").c_str());
    if (ret != 0)
        LOG << "FabberRunData::uname failed - uname.txt not created" << endl;

    if (ltl)
    {
        // try to make a link to the latest version. If this fails, it doesn't really matter.
        int ret = system(("ln -sfn '" + m_outdir + "' '" + basename + "_latest'").c_str());
        if (ret != 0)
            LOG << "FabberRunData::link-to-latest failed" << endl;
    }
#endif

    return m_outdir;
}

ostream &operator<<(ostream &out, const FabberRunData &opts)
{
    for (map<string, string>::const_iterator i = opts.m_params.begin(); i != opts.m_params.end();
         ++i)
    {
        if (i->second == "")
            out << "--" << i->first << endl;
        else
            out << "--" << i->first << "='" << i->second << "'" << endl;
    }
    return out;
}

const NEWMAT::Matrix &FabberRunData::GetMainVoxelData()
{
    // Main voxel data is a bit special because it can
    // come from multiple files
    try
    {
        return GetVoxelData("data");
    }
    catch (DataNotFound &e)
    {
        // See if we seem to have multi-data
        try
        {
            GetVoxelData("data1");
        }
        catch (DataNotFound &e2)
        {
            // Throw original exception
            throw(e);
        }
        return GetMainVoxelDataMultiple();
    }
}

const NEWMAT::Matrix &FabberRunData::GetVoxelSuppData()
{
    // FIXME Currently Fabber models assume that suppdata will return an empty matrix if none
    // present.
    try
    {
        return GetVoxelData("suppdata");
    }
    catch (DataNotFound &e)
    {
        return m_empty;
    }
}

int FabberRunData::GetVoxelDataSize(const std::string &key)
{
    NEWMAT::Matrix mat = GetVoxelData(key);
    return mat.Nrows();
}

const NEWMAT::Matrix &FabberRunData::GetVoxelCoords()
{
    return GetVoxelData("coords");
}

const NEWMAT::Matrix &FabberRunData::GetVoxelData(const std::string &key)
{
    // Attempt to load data if not already present. Will
    // throw an exception if parameter not specified
    // or file could not be loaded
    //
    // FIXME different exceptions? What about use case where
    // data is optional?
    string key_cur = key;
    string data_key = "";
    while (key_cur != "")
    {
        data_key = key_cur;
        key_cur = GetStringDefault(key_cur, "");
        // Avoid possible circular reference!
        if (key_cur == key)
            break;
    }

    const NEWMAT::Matrix &m = LoadVoxelData(data_key);
    return m;
}

const NEWMAT::Matrix &FabberRunData::LoadVoxelData(const std::string &key)
{
    if (m_voxel_data.count(key) == 0)
    {
        throw DataNotFound(key);
    }
    else
    {
        return m_voxel_data.find(key)->second;
    }
}

const Matrix &FabberRunData::GetMainVoxelDataMultiple()
{
    vector<Matrix> dataSets;
    int n = 1;
    while (true)
    {
        try
        {
            dataSets.push_back(GetVoxelData("data" + stringify(n)));
            n++;
        }
        catch (DataNotFound &e)
        {
            // No more data sets to combine, carry on with what we've got
            break;
        }
    }

    string order = GetStringDefault("data-order", "interleave");
    int nSets = dataSets.size();
    if (nSets < 1)
    {
        throw DataNotFound("data");
    }
    if ((order == "singlefile") && nSets > 1)
    {
        throw InvalidOptionValue("data-order", "singlefile", "More than one file specified");
    }

    if (order == "interleave")
    {
        LOG << "FabberRunData::Combining data into one big matrix by interleaving..." << endl;
        // Interleave - For example if the data sets are A, B, C and each
        // has 3 time points 1, 2, 3 the final time series will be
        // A1B1C1A2B2C2A3B3C3
        int nTimes = dataSets[0].Nrows();
        m_mainDataMultiple.ReSize(nTimes * nSets, dataSets[0].Ncols());
        for (int i = 0; i < nTimes; i++)
        {
            for (int j = 0; j < nSets; j++)
            {
                if (dataSets[j].Nrows() != nTimes)
                {
                    // Data sets need same number of time points if they are to be interleaved
                    throw InvalidOptionValue("data-order", "interleave",
                        "Data sets must all have the same number of time points");
                }
                m_mainDataMultiple.Row(nSets * i + j + 1) = dataSets.at(j).Row(i + 1);
            }
        }
    }
    else if (order == "concatenate")
    {
        LOG << "FabberRunData::Combining data into one big matrix by concatenating..." << endl;
        // Concatentate - For example if the data sets are A, B, C and each
        // has 3 time points 1, 2, 3 the final time series will be
        // A1A2A3B1B2B3C1C2C3
        m_mainDataMultiple = dataSets.at(0);
        for (unsigned j = 1; j < dataSets.size(); j++)
        {
            m_mainDataMultiple &= dataSets.at(j);
        }
    }
    else if (order == "singlefile")
    {
        m_mainDataMultiple = dataSets[0];
    }
    else
    {
        throw InvalidOptionValue("data-order", order, "Value not recognized");
    }

    LOG << "FabberRunData::Done loading data, size = " << m_mainDataMultiple.Nrows()
        << " timepoints by " << m_mainDataMultiple.Ncols() << " voxels" << endl;
    return m_mainDataMultiple;
}

void FabberRunData::ClearVoxelData(string key)
{
    if (key != "")
    {
        m_voxel_data.erase(key);
    }
    else
    {
        m_voxel_data.clear();
    }
}

void FabberRunData::SetVoxelData(string key, const NEWMAT::Matrix &data)
{
    CheckSize(key, data);
    m_voxel_data[key] = data;
}

void FabberRunData::SaveVoxelData(
    const std::string &filename, NEWMAT::Matrix &data, VoxelDataType data_type)
{
    LOG << "FabberRunData::Saving to memory: " << filename << endl;
    // FIXME what should we do with data_type?
    SetVoxelData(filename, data);
}

void FabberRunData::SetVoxelCoords(const NEWMAT::Matrix &coords)
{
    // We assume 3D coordinates. Fabber could work for different
    // numbers of dimensions but would require extensive refactoring
    if (coords.Ncols() > 0 && coords.Nrows() != 3)
    {
        throw InvalidOptionValue("Coordinates dimensions", stringify(coords.Nrows()),
            "Co-ordinates must be 3 dimensional");
    }

    SetVoxelData("coords", coords);

    if (m_extent.size() == 0)
    {
        m_extent.resize(3);
        m_dims.resize(3);
        for (int i = 0; i < 3; i++)
        {
            if (coords.Ncols() > 0)
            {
                // FIXME we assume coords will not be negative
                m_extent[i] = coords.Row(i + 1).Maximum() - coords.Row(i + 1).Minimum() + 1;
            }
            else
            {
                m_extent[i] = 0;
            }
            m_dims[i] = 1.0;
        }
    }
}

void FabberRunData::GetExtent(vector<int> &extent, vector<float> &dims)
{
    extent = m_extent;
    dims = m_dims;
}

void FabberRunData::SetExtent(int nx, int ny, int nz, float sx, float sy, float sz)
{
    if ((nx < 0) || (ny < 0) || (nz < 0) || (sx <= 0) || (sy <= 0) || (sz <= 0))
    {
        throw InvalidOptionValue(
            "extent", "negative values", "Extent and voxel sizes must be positive");
    }

    m_extent.resize(3);
    m_dims.resize(3);
    m_extent[0] = nx;
    m_extent[1] = ny;
    m_extent[2] = nz;
    m_dims[0] = sx;
    m_dims[1] = sy;
    m_dims[2] = sz;
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

vector<vector<int> > &FabberRunData::GetNeighbours(int n_dims)
{
    if (m_neighbours.size() > 0)
        return m_neighbours;

    const Matrix &coords = GetVoxelCoords();
    const int nvoxels = coords.Ncols();
    if (nvoxels == 0)
        return m_neighbours;

    // Voxels must be ordered by increasing z, y and x values respectively
    // otherwise binary search for voxel by offset will not work
    // CheckCoordMatrixCorrectlyOrdered(coords);

    // Create a column vector with one entry per voxel.
    ColumnVector offsets(nvoxels);

    // Populate offsets with the offset into the
    // matrix of each voxel. We assume that co-ordinates
    // could be zero but not negative
    int xsize = coords.Row(1).Maximum() + 1;
    int ysize = coords.Row(2).Maximum() + 1;
    for (int v = 1; v <= nvoxels; v++)
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
    // For example if n_dims=2, max_delta=3 so we
    // only look for neighbours in rows and columns
    //
    // However note we still need the full list of 3D deltas for later
    int max_delta = n_dims * 2 - 1;

    // Neighbours is a vector of vectors, so each voxel
    // will have an entry which is a vector of its neighbours
    m_neighbours.resize(nvoxels);

    // Go through each voxel. Note that offsets is indexed from 1 not 0
    // however the offsets themselves (potentially) start at 0.
    for (int vid = 1; vid <= nvoxels; vid++)
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
            // Changed (fixed)? because if n_dims != 3 we still need
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
            m_neighbours.at(vid - 1).push_back(id);
        }
    }
    return m_neighbours;
}

vector<vector<int> > &FabberRunData::GetSecondNeighbours(int n_dims)
{
    if (m_neighbours2.size() > 0)
        return m_neighbours2;

    GetNeighbours(n_dims);
    const int nvoxels = m_neighbours.size();
    if (nvoxels == 0)
        return m_neighbours2;

    // Similar algorithm but looking for Neighbours-of-neighbours, excluding self,
    // but including duplicates if there are two routes to get there
    // (diagonally connected)
    m_neighbours2.resize(nvoxels);
    for (int vid = 1; vid <= nvoxels; vid++)
    {
        // Go through the list of neighbours for each voxel.
        for (unsigned n1 = 0; n1 < m_neighbours.at(vid - 1).size(); n1++)
        {
            // n1id is the voxel index (not the offset) of the neighbour
            int n1id = m_neighbours[vid - 1].at(n1);
            int checkNofN = 0;
            // Go through each of it's neighbours. Add each, apart from original voxel
            for (unsigned n2 = 0; n2 < m_neighbours.at(n1id - 1).size(); n2++)
            {
                int n2id = m_neighbours[n1id - 1].at(n2);
                if (n2id != vid)
                {
                    m_neighbours2[vid - 1].push_back(n2id);
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

    return m_neighbours2;
}

void FabberRunData::CheckSize(std::string key, const NEWMAT::Matrix &mat)
{
    if (m_voxel_data.size() > 0)
    {
        int nvoxels = m_voxel_data.begin()->second.Ncols();
        if (mat.Ncols() != nvoxels)
        {
            throw InvalidOptionValue("Voxels in " + key, stringify(mat.Ncols()),
                "Incorrect size - should contain " + stringify(nvoxels));
        }
    }
}
