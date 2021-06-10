// fwdmodel_exp.h - A simple exponential decay fitting model
#pragma once

#include "fabber_core/fwdmodel.h"

#include "armawrap/newmat.h"

#include <string>
#include <vector>

class ExpFwdModel : public FwdModel {
public:
    static FwdModel* NewInstance();

    ExpFwdModel()
        : m_num(1), m_dt(1.0)
    {
    }

    std::string ModelVersion() const;
    std::string GetDescription() const;
    void GetOptions(std::vector<OptionSpec> &opts) const;

    void Initialize(FabberRunData &args);
    void EvaluateModel(const NEWMAT::ColumnVector &params,
                       NEWMAT::ColumnVector &result,
                       const std::string &key="") const;

    void InitVoxelPosterior(MVNDist &posterior) const;

protected:
    void GetParameterDefaults(std::vector<Parameter> &params) const;

private:
    int m_num;
    double m_dt;
    static FactoryRegistration<FwdModelFactory, ExpFwdModel> registration;
};
