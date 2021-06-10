/*  noisemodel_ar.h - Class declaration for the multiple white noise model

 Adrian Groves and Michael Chappell, FMRIB Image Analysis Group & IBME QuBIc Group

 Copyright (C) 2007-2015 University of Oxford  */

/*  CCOPYRIGHT */
#pragma once

#include "noisemodel.h"

#include "dist_gamma.h"

#include "armawrap/newmat.h"

#include <ostream>
#include <string>
#include <vector>

/**
 * Parameters for the white noise model
 *
 * The noise model can have any number of parameters (Phis)
 * which apply to different samples in the timeseries, e.g.
 * a noise pattern of 12121212... has two parameters
 * one applying to odd numbered samples, and one to even.
 *
 * Each Phi is associated with a gamma distribution
 */
class WhiteParams : public NoiseParams
{
public:
    explicit WhiteParams(int N);
    WhiteParams(const WhiteParams &from);
    virtual const WhiteParams &operator=(const NoiseParams &in);

    virtual WhiteParams *Clone() const;

    virtual const MVNDist OutputAsMVN() const;
    virtual void InputFromMVN(const MVNDist &mvn);
    virtual void Dump(std::ostream &os) const;

private:
    friend class WhiteNoiseModel;
    const int nPhis;
    std::vector<GammaDist> phis;
};

class WhiteNoiseModel : public NoiseModel
{
public:
    /**
     * Create a new instance of white noise. Used by factory
     * to create noise models by name. Parameters are set
     * during initialization
     */
    static NoiseModel *NewInstance();

    virtual void Initialize(FabberRunData &args);
    virtual WhiteParams *NewParams() const;
    int NumParams();

    virtual void HardcodedInitialDists(NoiseParams &prior, NoiseParams &posterior) const;

    virtual void UpdateNoise(NoiseParams &noise, const NoiseParams &noisePrior,
        const MVNDist &theta, const LinearFwdModel &model, const NEWMAT::ColumnVector &data) const;

    virtual void UpdateTheta(const NoiseParams &noise, MVNDist &theta, const MVNDist &thetaPrior,
        const LinearFwdModel &model, const NEWMAT::ColumnVector &data,
        MVNDist *thetaWithoutPrior = NULL, float LMalpha = 0) const;

    virtual double CalcFreeEnergy(const NoiseParams &noise, const NoiseParams &noisePrior,
        const MVNDist &theta, const MVNDist &thetaPrior, const LinearFwdModel &model,
        const NEWMAT::ColumnVector &data) const;

protected:
    /** Pattern of noise distributions as they apply to points in time series */
    std::string phiPattern;

    /** Allow phi to be locked externally */
    double lockedNoiseStdev;

    /** Allow external setting of the prior noise std deviation (and thence phi) */
    double phiprior;

    /**
     * Diagonal matrices, indicating which data points use each phi
     *
     * Mutable because it's initialized lazily by MakeQis
     */
    mutable std::vector<NEWMAT::DiagonalMatrix> Qis;

    /** Create Qis vector */
    void MakeQis(int dataLen) const;
};
