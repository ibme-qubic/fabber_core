/*  fwdmodel_poly.h - Implements a simple polynomial model

 Copyright (C) 2007 University of Oxford  */

/*  CCOPYRIGHT */
#pragma once

#include "dist_mvn.h"
#include "fwdmodel.h"
#include "rundata.h"

#include <newmat.h>

#include <string>
#include <vector>

/**
 * Forward model which fits to a simple polynomial function
 * e.g. a + bx + cx^2 + ....
 *
 * The parameter 'degree' can be used to set the maximum power
 * in the polynomial, e.g. if degree=3, there will be 4 parameters
 *
 * Note that this class is mostly for testing purposes and is not
 * designed to be overriden.
 */
class PolynomialFwdModel : public FwdModel
{
public:
    static FwdModel *NewInstance();

    PolynomialFwdModel()
        : m_degree(0)
    {
    }
    void GetOptions(std::vector<OptionSpec> &opts) const;
    std::string GetDescription() const;
    std::string ModelVersion() const;

    void Initialize(FabberRunData &args);
    int NumParams() const;
    void NameParams(std::vector<std::string> &names) const;

    void HardcodedInitialDists(MVNDist &prior, MVNDist &posterior) const;
    void EvaluateModel(const NEWMAT::ColumnVector &params, NEWMAT::ColumnVector &result,
        const std::string &key = "") const;

private:
    int m_degree;
};
