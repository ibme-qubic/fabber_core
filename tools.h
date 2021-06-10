/*  tools.h - miscellaneous useful functions & classes

 Adrian Groves, FMRIB Image Analysis Group

 Copyright (C) 2007 University of Oxford  */

/*  CCOPYRIGHT */
#pragma once

#include "rundata.h"

#include "armawrap/newmat.h"

#include <assert.h>
#include <math.h>

namespace fabber
{
/**
 * Read 'small' matrix from file.
 *
 * The matrix may be in ASCII or VEST format
 */
NEWMAT::Matrix read_matrix_file(std::string filename);

NEWMAT::ReturnMatrix MaskRows(NEWMAT::Matrix m, std::vector<int> masked_rows);
NEWMAT::ReturnMatrix MaskRows(NEWMAT::ColumnVector v, std::vector<int> masked_rows);
}

// Calculate log-gamma from a Taylor expansion; good to one part in 2e-10.
double gammaln(double x);

/**
 * Base class for a generic 1-dimensional function which takes a double and returns a double
 *
 * Could add more functions.. e.g. initialGuess, domain of validity, etc.
 */
class GenericFunction1D : public Loggable
{
public:
    /**
     * Calculate the value of this function
     *
     * @param x Input value
     * @return output of function
     */
    virtual double Calculate(double x) const = 0;

    /**
     * Allow us to calculate the value of the function
     * using the () operator
     *
     * e.g. MyFunction myfunc;
     *      double answer = myfunc(4);
     */
    double operator()(double x) const
    {
        return Calculate(x);
    }
    /**
     * This is useful if your function is very slow to calculate, but you have
     * some cached partial calculations available.  If you have a suitable
     * cached value, store it into guess and return true.  Otherwise return
     * false (and leave guess unchanged).
     */
    virtual bool PickFasterGuess(
        double *guess, double lower, double upper, bool allowEndpoints = false) const
    {
        return false;
    }

    virtual ~GenericFunction1D()
    {
    }

private:
    // Function's constant data should go here
};

#define REALMAX (1.7976931348623158e+308)

/**
 * Base class for a method of guessing the next value to try if we are looking
 * for the zero of a function
 *
 * @param lower Lower value of input
 * @param upper Upper value of input. Must have lower<upper
 * @param atLower Value of function at lower input
 * @param atUpper Value of function at upper input
 * @return Next value to try
 */
class Guesstimator : public Loggable
{
public:
    virtual double GetGuess(double lower, double upper, double atLower, double atUpper) = 0;
    virtual ~Guesstimator()
    {
    }
};

/**
 * Guesstimator which just suggests you try half way in between the two values
 */
class BisectionGuesstimator : public Guesstimator
{
public:
    virtual double GetGuess(double lower, double upper, double, double)
    {
        assert(lower < upper);
        return (lower + upper) / 2;
    }
};

/**
 * Guesstimator which suggests the geometric mean of lower and upper
 */
class LogBisectionGuesstimator : public Guesstimator
{
public:
    virtual double GetGuess(double lower, double upper, double, double)
    {
        assert(lower > 0 && upper > lower);
        double guess = sqrt(lower * upper);
        if (lower >= guess || guess >= upper)
        {
            cout << "Uh-oh... lower = " << lower << ", guess = " << guess << ", upper = " << upper
                 << endl;
        }
        return guess;
    }
};

/**
 * equations below: from NRIC, section 9.2.  Simpler than Brent, slightly less reliable.
 */
class RiddlersGuesstimator : public Guesstimator
{
public:
    virtual double GetGuess(double lower, double upper, double atLower, double atUpper);
    RiddlersGuesstimator()
        : halfDone(false)
        , x1(0)
        , x2(0)
        , fx1(0)
        , fx2(0)
    {
    }

private:
    bool halfDone;           // waiting for f(x3) result?
    double x1, x2, fx1, fx2; // save from phase 1 for phase 2; only valid when halfDone is true.
};

/**
 * Performs a log transformation on a RiddlesGuesstimator
 */
class LogRiddlersGuesstimator : public RiddlersGuesstimator
{
public:
    virtual double GetGuess(double lower, double upper, double atLower, double atUpper)
    {
        return exp(RiddlersGuesstimator::GetGuess(log(lower), log(upper), atLower, atUpper));
    }
};

/**
 * Performs linear interpolation to produce the next guess
 */
class InterpGuesstimator : public Guesstimator
{
public:
    virtual double GetGuess(double lower, double upper, double atLower, double atUpper)
    {
        double guess = upper - atUpper * (upper - lower) / (atUpper - atLower);

        if (lower >= guess || guess >= upper)
        {
            cout << "Uh-oh... lower = " << lower << ", guess = " << guess << ", upper = " << upper
                 << ", atLower = " << atLower << ", atUpper = " << atUpper << endl;
        }
        return guess;
    }
};

/**
 * Finds the zero of a GenericFunction
 *
 * Note that you have to specify one or more tolerances to get any sensible
 * results.  Also note that the ratio tolerances current assume that the
 * X or Y value always positive -- otherwise it'll stop too early!
 */
class ZeroFinder : public Loggable
{
public:
    explicit ZeroFinder(const GenericFunction1D &f)
        : fcn(f)
        , searchMin(-REALMAX)
        , searchMax(REALMAX)
        , searchGuess(0)
        , searchScale(REALMAX)
        , searchScaleGrowth(2)
        , maxEvaluations(1000000)
        , tolX(REALMAX)
        , tolY(REALMAX)
        , ratioTolX(REALMAX)
        , ratioTolY(REALMAX)
        , guesstimator(new BisectionGuesstimator())
        , verbosity(2)
    {
        m_log = f.GetLogger();
    }
    /**
     * Return input value at which function is zero
     */
    virtual double FindZero() const = 0;

    /**
     * Returns the input value at which the function is zero,
     * using the () operator
     *
     * e.g. ZeroFinder finder(MyFunc);
     *      root = finder();
     */
    operator double() const
    {
        return FindZero();
    }
    virtual ~ZeroFinder()
    {
    }
    /**
     * Set initial guess
     */
    ZeroFinder &InitialGuess(double guess)
    {
        searchGuess = guess;
        return *this;
    }
    /**
     * Set the a minimum value we will not search below
     */
    ZeroFinder &SearchMin(double min)
    {
        searchMin = min;
        return *this;
    }
    /**
     * Set the a maximum value we will not search above
     */
    ZeroFinder &SearchMax(double max)
    {
        searchMax = max;
        return *this;
    }
    ZeroFinder &InitialScale(double scale)
    {
        searchScale = scale;
        return *this;
    }
    ZeroFinder &ScaleGrowth(double growth)
    {
        assert(growth > 1);
        searchScaleGrowth = growth;
        return *this;
    }
    /**
     * Set the a maximum number of trials before we stop
     */
    ZeroFinder &MaxEvaluations(int evals)
    {
        assert(evals > 1);
        maxEvaluations = evals;
        return *this;
    }
    ZeroFinder &TolX(double tx)
    {
        assert(tx > 0);
        tolX = tx;
        return *this;
    }
    ZeroFinder &TolY(double ty)
    {
        assert(ty > 0);
        tolY = ty;
        return *this;
    }
    //    ZeroFinder& RatioTolY(double rty) // utterly pointless -- looking for a sign change!
    //      { assert(rty>1); ratioTolY = rty; return *this; }
    ZeroFinder &RatioTolX(double rtx)
    {
        assert(rtx > 1);
        ratioTolX = rtx;
        return *this;
    }
    /**
     * Set a Guesstimator to use to produce the next estimate
     */
    ZeroFinder &SetGuesstimator(Guesstimator *g)
    {
        delete guesstimator;
        guesstimator = g;
        return *this;
    }
    ZeroFinder &Verbosity(int v)
    {
        verbosity = v;
        return *this;
    }

protected:
    const GenericFunction1D &fcn;

    // Optional parameters:
    double searchMin;
    double searchMax;
    double searchGuess;
    double searchScale;
    double searchScaleGrowth;
    int maxEvaluations;
    double tolX;
    double tolY;
    double ratioTolX;
    double ratioTolY;
    Guesstimator *guesstimator;
    int verbosity;
};

class DescendingZeroFinder : public ZeroFinder
{
public:
    explicit DescendingZeroFinder(const GenericFunction1D &f)
        : ZeroFinder(f)
    {
        return;
    }
    virtual double FindZero() const;
};
