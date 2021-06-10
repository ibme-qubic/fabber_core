/**
 * prior.h
 *
 * Classes for parameter priors
 *
 * Copyright (C) 2007-2017 University of Oxford
 */

#include "dist_mvn.h"
#include "fwdmodel.h"
#include "run_context.h"
#include "rundata.h"

#include "armawrap/newmat.h"

#include <ostream>
#include <string>
#include <vector>

const char PRIOR_NORMAL = 'N';    // non-spatial prior
const char PRIOR_IMAGE = 'I';     // image prior
const char PRIOR_ARD = 'A';       // ARD prior
const char PRIOR_SPATIAL_M = 'M'; // Markov random field - normally used
const char PRIOR_SPATIAL_m = 'm'; // 'M' with Dirichlet BCs
const char PRIOR_SPATIAL_P = 'P'; // Alternative to M (Penny prior?)
const char PRIOR_SPATIAL_p = 'p'; // P with Dirichlet BCs
const char PRIOR_DEFAULT = '-';   // Use whatever the model specifies

/**
 * Abstract interface for a parameter prior
 */
class Prior : public Loggable
{
public:
    virtual ~Prior()
    {
    }
    /** Dump info to output stream */
    virtual void DumpInfo(std::ostream &out) const = 0;

    /**
     * Apply prior information to an MVN
     *
     * Returns any additional free energy contribution (e.g. for ARD priors)
     */
    virtual double ApplyToMVN(MVNDist *prior, const RunContext &ctx) = 0;

    /** Expand the param-spatial-priors string to give a value for each parameter */
    static std::string ExpandPriorTypesString(std::string priors_str, unsigned int num_params);
};

/**
 * Prior which has a mean and precision
 */
class DefaultPrior : public Prior
{
public:
    DefaultPrior(const Parameter &param);
    virtual ~DefaultPrior()
    {
    }
    /** Parameter name this prior applies to */
    std::string m_param_name;

    /** Parameter index number */
    unsigned int m_idx;

    /** Prior type code */
    char m_type_code;

    /** Prior mean and variance*/
    DistParams m_params;

    virtual void DumpInfo(std::ostream &out) const;
    virtual double ApplyToMVN(MVNDist *prior, const RunContext &ctx);
};

/**
 * Prior which has takes its mean from a per-voxel image, with a constant precision
 */
class ImagePrior : public DefaultPrior
{
public:
    ImagePrior(const Parameter &param, FabberRunData &rundata);

    virtual void DumpInfo(std::ostream &out) const;
    virtual double ApplyToMVN(MVNDist *prior, const RunContext &ctx);

protected:
    /** Filename containing image data if required */
    std::string m_filename;

    /** Image data if required */
    NEWMAT::RowVector m_image;
};

/**
 * ARD prior
 */
class ARDPrior : public DefaultPrior
{
public:
    ARDPrior(const Parameter &param, FabberRunData &rundata)
        : DefaultPrior(param)
    {
        m_log = rundata.GetLogger();
    }

    virtual void DumpInfo(std::ostream &out) const;
    virtual double ApplyToMVN(MVNDist *prior, const RunContext &ctx);
};

/**
 * Prior which uses spatial information to inform the prior
 *
 * Currently this merges all different spatial priors into one class. Ideally
 * find a way to split each one into a separate subclass.
 */
class SpatialPrior : public DefaultPrior
{
public:
    SpatialPrior(const Parameter &param, FabberRunData &rundata);

    virtual void DumpInfo(std::ostream &out) const;
    virtual double ApplyToMVN(MVNDist *prior, const RunContext &ctx);

protected:
    double CalculateaK(const RunContext &ctx);
    double m_aK;
    int m_spatial_dims;
    double m_spatial_speed;
    bool m_update_first_iter;
    double m_q1;
    double m_q2;
};

/**
 * Creates instances of Prior depending on the input options
 */
class PriorFactory : public Loggable
{
public:
    PriorFactory(FabberRunData &rundata);

    /** Create priors for all model parameters */
    std::vector<Prior *> CreatePriors(const std::vector<Parameter> &params);

private:
    FabberRunData &m_rundata;

    /** Create a prior for a parameter */
    Prior *CreatePrior(Parameter p);
};

std::ostream &operator<<(std::ostream &out, const Prior &value);
