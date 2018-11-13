/*  inference.h - General inference technique base class

 Adrian Groves & Michael Chappell, FMRIB Image Analysis Group & IBME QuBIc Group

 Copyright (C) 2007-2015 University of Oxford  */

/*  CCOPYRIGHT */

#pragma once

#include "dist_mvn.h"
#include "easylog.h"
#include "factories.h"
#include "fwdmodel.h"
#include "noisemodel.h"
#include "rundata.h"
#include "run_context.h"

#include <iomanip>
#include <string>
#include <vector>

class InferenceTechnique : public Loggable
{
public:
    /**
     * Static member function to return the names of all known
     * inference techniques
     */
    static std::vector<std::string> GetKnown();

    /**
     * Static member function, to pick an inference technique from a name
     */
    static InferenceTechnique *NewFromName(const std::string &name);

    /**
     * Get usage information for a named method
     */
    static void UsageFromName(const std::string &name, std::ostream &stream);

    /**
     * Default constructor.
     */
    InferenceTechnique();

    /**
     * Get option descriptions for this inference method.
     */
    virtual void GetOptions(std::vector<OptionSpec> &opts) const {};

    /**
     * @return human-readable description of the inference method.
     */
    virtual std::string GetDescription() const = 0;

    /**
     * Get the code version. There is no fixed format for this,
     * and it has no meaning other than by comparison with different
     * versions of the same inference method code.
     *
     * See fwdmodel.cc for an example of how to implement this to
     * return a CVS file version.
     *
     * @return a string indicating the inference code version.
     */
    virtual std::string GetVersion() const = 0;

    /**
     * Initialize a new instance to use the given forward model
     * and extract additional configuration from the given
     * arguments.
     * @param fwd_model Forward model to be used.
     * @param args Additional configuration parameters.
     */
    virtual void Initialize(FabberRunData &rundata);

    /**
     * Perform inference using the given model upon the given data.
     *
     * This method should only be called after Initialize()
     * Subclasses of InferenceTechnique must implement this method
     * to carry out their given inference calculations
     */
    virtual void DoCalculations(FabberRunData &rundata) = 0;

    /**
     * Save the results
     */
    virtual void SaveResults(FabberRunData &rundata) const;

    /**
     * Destructor.
     */
    virtual ~InferenceTechnique();

protected:
    /** 
     * Stores current run state (parameters, MVNs, linearization centres etc 
     */
    std::vector<ThreadContext *> m_ctxs;

    /**
     * If true, stop if we get a numerical execption in any voxel. If false,
     * simply print a warning and continue
     */
    bool m_halt_bad_voxel;

    /**
     * Include very verbose debugging output
     */
    bool m_debug;

private:
    /**
     * Private to prevent assignment
     */
    const InferenceTechnique &operator=(const InferenceTechnique &from)
    {
        assert(false);
        return from;
    }
};

/**
 * \ref SingletonFactory that returns pointers to
 * \ref InferenceTechnique.
 */
typedef SingletonFactory<InferenceTechnique> InferenceTechniqueFactory;
