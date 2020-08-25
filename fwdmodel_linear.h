/*  fwdmodel_linear.h - Linear forward model and related classes

 Adrian Groves & Michael Chappell, FMRIB Image Analysis Group & IBME QuBIc Group

 Copyright (C) 2007-20015 University of Oxford  */

/*  CCOPYRIGHT */

#pragma once

#include "fwdmodel.h"

#include "armawrap/newmat.h"

#include <string>
#include <vector>

/**
 * Linear forward model
 *
 * This is a model in which the output depends linearly on
 * the parameters
 *
 * The model is defined by a design matrix whose columns are each
 * a timeseries for a single parameter. The number of rows must
 * therefore equal the number of timeseries samples.
 */
class LinearFwdModel : public FwdModel
{
public:
    static FwdModel *NewInstance();

    virtual void GetOptions(std::vector<OptionSpec> &opts) const;
    virtual std::string GetDescription() const;
    virtual std::string ModelVersion() const;

    virtual void Initialize(FabberRunData &args);

    /**
     * Evaluate the model.
     *
     * In the simplest case the model is evaluated by taking the product of the
     * design
     * matrix with the parameter vector. However the model can be recentred and
     * offsetted
     * so with this in mind the evaluation consists of:
     *
     * R = J x (P - C) + O
     *
     * Where R is the output result, J is the Jacobian or design matrix, P are the
     * parameters,
     * C is the centre and O is the offset.
     *
     * Centre and offset will be zero for a basic linear model, however a subclass
     * can set
     * these to non-zero values, e.g. in LinearizedFwdModel
     *
     * @param params Parameter vector, P in above equation. Length equal to number
     * of parameters
     *               i.e. number of columns in design matrix.
     * @param result Result vector, R in above equation. Length equal to number of
     * time samples
     *               i.e. number of rows in design matrix.
     */
    virtual void EvaluateModel(const NEWMAT::ColumnVector &params, NEWMAT::ColumnVector &result,
        const std::string &key = "") const;

    /**
     * @return the Jacobian, or design matrix
     */
    NEWMAT::ReturnMatrix Jacobian() const;

    /**
     * @return the vector used to recentre the parameters
     */
    NEWMAT::ReturnMatrix Centre() const;

    /**
     * @return the vector used to offset the result vector
     */
    NEWMAT::ReturnMatrix Offset() const;

protected:
    virtual void GetParameterDefaults(std::vector<Parameter> &params) const;

    NEWMAT::Matrix m_jacobian;     // J (tranposed?)
    NEWMAT::ColumnVector m_centre; // m
    NEWMAT::ColumnVector m_offset; // g(m) - The amount to effectively subtract from Y is g(m)-J*m
};

/**
 * Linearized wrapper interface to another nonlinear forward model
 *
 * This is not used as a model directly so we do not need to implement
 * parameter defaults, version numbers, options, etc. It is debatable whether
 * it should be a FwdModel instance at all.
 */
class LinearizedFwdModel : public LinearFwdModel
{
public:
    /**
     * Create a linear approximation to another (nonlinear) forward model
     *
     * This model implements a linear approximation to the given
     * model. The ReCentre method must be called to update the
     * Jacobian, offset and centre so that Evaluate returns
     * the linear approximation which is correct at the given
     * centre.
     *
     * The model pointer is not owned by this class and will not be freed
     */
    explicit LinearizedFwdModel(const FwdModel *model);

    /**
     * Copy constructor (needed for using vector<LinearizedFwdModel>)
     *
     * NOTE: This is a reference, not a pointer... and it *copies* the
     * given LinearizedFwdModel, rather than using it as its nonlinear model!
     */
    LinearizedFwdModel(const LinearizedFwdModel &from);

    /**
     * Re-calculate the linearized model about the given centre
     *
     * @param about New centre in parameter space
     *
     * This function sets the model centre to the vector given
     * and recalculates the Jacobian and offset to form a linear
     * approximation to the underlying model.
     *
     * The offset is set to the evaluation of the underlying
     * model on the new centre, to ensure that Evaluate is
     * correct when called with the given centre
     *
     * The Jacobian is calculated either by asking the model
     * directly using the \ref Gradient method, or if that is not
     * implemented by numerical differentiation about the
     * new centre
     */
    void ReCentre(const NEWMAT::ColumnVector &about);

private:
    const FwdModel *m_model;
};
