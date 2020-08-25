/*  noisemodel.h - Class declaration for generic noise models

 Adrian Groves & Michael Chappell, FMRIB Image Analysis Group & IBME QuBIc Group

 Copyright (C) 2007-2015 University of Oxford  */

/*  CCOPYRIGHT */

#pragma once

#include "dist_mvn.h"
#include "factories.h"
#include "fwdmodel_linear.h"
#include "rundata.h"

#include "armawrap/newmat.h"

#include <ostream>
#include <string>

/**
 * Base class for parameters to a noise model
 *
 * Each derived NoiseModel will have a derived NoiseParams
 */
class NoiseParams
{
public:
    /** Virtual destructor for subclasses */
    virtual ~NoiseParams()
    {
    }
    /**
     * Create a copy of the params
     *
     * Polymorphic so subclasses can create copies of the correct class
     */
    virtual NoiseParams *Clone() const = 0;

    /** Assignment operator, needs to check subclass equivalence */
    virtual const NoiseParams &operator=(const NoiseParams &in) = 0;

    /**
     * Output as a multivariate normal dist. Noise models need not
     * be normal, but they will generally have means and precisions
     */
    virtual const MVNDist OutputAsMVN() const = 0;

    /**
     * Input from MVN distribution - see OutputAsMVN
     */
    virtual void InputFromMVN(const MVNDist &mvn) = 0;

    /**
     * Dump human-readable debug output to output stream
     */
    virtual void Dump(std::ostream &os) const = 0;
};

/**
 * Generic noise model
 *
 * This class & derived classes should be essentially data-free, instead storing
 * the relevant noise parameters in a NoiseParams-derived subclass
 *
 * FIXME unclear why params need to be separate?
 */
class NoiseModel : public Loggable
{
public:
    /**
     * Static member function, to pick a noise model from a name
     */
    static NoiseModel *NewFromName(const std::string &name);

    /**
     * Create a new instance of this class. Subclasses
     * implement this to produce an instance of themselves
     *
     * @return pointer to new instance.
     */
    static NoiseModel *NewInstance();

    /** Virtual destrictor for subclasses */
    virtual ~NoiseModel()
    {
    }
    /**
     * Initialize a new instance using configuration from the given
     * arguments.
     *
     * @param args Configuration parameters.
     */
    virtual void Initialize(FabberRunData &args);

    /**
     * Create a new instance of the noise parameters
     */
    virtual NoiseParams *NewParams() const = 0;

    /**
     * Suggest some nice default values for noise parameters:
     */
    virtual void HardcodedInitialDists(NoiseParams &prior, NoiseParams &posterior) const = 0;

    /**
     * Some noise models might want to precalculate things (for efficiency
     * reasons), based on the length of the data... if you don't know what
     * this is for then just ignore it.
     */
    virtual void Precalculate(NoiseParams &noise, const NoiseParams &noisePrior,
        const NEWMAT::ColumnVector &sampleData) const
    {
    }

    // VB Updates
    // The following could potentially be split into substeps; but since
    // these would necessarily be model-specific, it's nice to have a
    // general catch-all update step.  Presumably this function
    // would call all the other functions in some order.

    /**
     * Update noise parameters
     */
    virtual void UpdateNoise(NoiseParams &noise, const NoiseParams &noisePrior,
        const MVNDist &theta, const LinearFwdModel &model,
        const NEWMAT::ColumnVector &data) const = 0;

    /**
     * Update model parameters
     *
     * @param thetaWithoutPrior used for --spatial-prior-output-correction
     */
    virtual void UpdateTheta(const NoiseParams &noise, MVNDist &theta, const MVNDist &thetaPrior,
        const LinearFwdModel &model, const NEWMAT::ColumnVector &data,
        MVNDist *thetaWithoutPrior = NULL, float LMalpha = 0) const = 0;

    /**
     * Calculate the free energy
     *
     * This may be used to assess convergence
     */
    virtual double CalcFreeEnergy(const NoiseParams &noise, const NoiseParams &noisePrior,
        const MVNDist &theta, const MVNDist &thetaPrior, const LinearFwdModel &model,
        const NEWMAT::ColumnVector &data) const = 0;

    /**
     * Initialize ARD
     */
    double SetupARD(vector<int> ardindices, const MVNDist &theta, MVNDist &thetaPrior) const;

    /**
     * Update ARD
     */
    double UpdateARD(vector<int> ardindices, const MVNDist &theta, MVNDist &thetaPrior) const;

    /**
     * Return the number of noise parameters
     */
    virtual int NumParams() = 0;

protected:
    /**
     * List of masked timepoints
     *
     * Masked timepoints are indexed starting at 1 and are ignored
     * in the analysis and parameter updates.
     */
    std::vector<int> m_masked_tpoints;

private:
    /**
     * Prevent copying using anything other than the Clone() function.
     * Could implement it, but not particularly useful and the default
     * shallow copy is not right.
     */
    const NoiseModel &operator=(const NoiseModel &) const
    {
        assert(false);
        return *this;
    }
};

/**
 * \ref SingletonFactory that returns pointers to \ref NoiseModel.
 */
typedef SingletonFactory<NoiseModel> NoiseModelFactory;
