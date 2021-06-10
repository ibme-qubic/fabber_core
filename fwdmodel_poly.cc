/*  fwdmodel_poly.cc - Implements a polynomial model

 Martin Craig

 Copyright (C) 2016 University of Oxford  */

/*  CCOPYRIGHT */

#include "fwdmodel_poly.h"

#include "dist_mvn.h"
#include "rundata.h"
#include "version.h"

#include "armawrap/newmat.h"

#include <string>
#include <vector>

using namespace std;

static OptionSpec OPTIONS[]
    = { { "degree", OPT_INT, "Maximum power in the polynomial function", OPT_REQ, "" }, { "" } };

FwdModel *PolynomialFwdModel::NewInstance()
{
    return new PolynomialFwdModel();
}
void PolynomialFwdModel::GetOptions(vector<OptionSpec> &opts) const
{
    for (int i = 0; OPTIONS[i].name != ""; i++)
    {
        opts.push_back(OPTIONS[i]);
    }
}

string PolynomialFwdModel::GetDescription() const
{
    return "Model which fits data to a simple polynomial function: c0 + c1x + "
           "c2x^2 ... etc";
}

string PolynomialFwdModel::ModelVersion() const
{
    return fabber_version();
}
void PolynomialFwdModel::Initialize(FabberRunData &args)
{
    FwdModel::Initialize(args);
    m_degree = convertTo<int>(args.GetString("degree"));
}

void PolynomialFwdModel::GetParameterDefaults(std::vector<Parameter> &params) const
{
    for (int i = 0; i < m_degree + 1; i++)
    {
        params.push_back(
            Parameter(i, "c" + stringify(i), DistParams(0, 1e12), DistParams(0, 1e12)));
    }
}

void PolynomialFwdModel::EvaluateModel(
    const NEWMAT::ColumnVector &params, NEWMAT::ColumnVector &result, const std::string &key) const
{
    assert(params.Nrows() == m_degree + 1);
    result.ReSize(data.Nrows());

    // LOG << "Evaluating for : " << params.t() << endl;
    for (int i = 1; i <= result.Nrows(); i++)
    {
        double res = 0;
        int p = 1;
        for (int n = 0; n <= m_degree; n++)
        {
            res += params(n + 1) * p;
            p *= i;
        }
        result(i) = res;
    }
}
