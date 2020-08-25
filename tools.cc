/*  tools.cc - miscellaneous useful functions & classes

 Adrian Groves, FMRIB Image Analysis Group

 Copyright (C) 2007 University of Oxford  */

/*  CCOPYRIGHT */

#include "tools.h"

#include "easylog.h"

#include <miscmaths/miscmaths.h>
#include "armawrap/newmat.h"

#include <limits>

using namespace std;
using NEWMAT::Matrix;
using NEWMAT::ColumnVector;
using MISCMATHS::read_vest;
using MISCMATHS::read_ascii_matrix;

namespace fabber
{
Matrix read_matrix_file(std::string filename)
{
    // Detect if file contains a VEST matrix or not
    try
    {
        return read_vest(filename);
    }
    catch (...)
    {
        // Do not care why this failed, if it was 'file not found'
        // we will discover that now we try to read it as ASCII
        return read_ascii_matrix(filename);
    }
}

ReturnMatrix MaskRows(Matrix m, vector<int> masked_rows)
{
    if (masked_rows.size() == 0)
    {
        return m;
    }
    else
    {
        NEWMAT::Matrix masked(m.Nrows() - masked_rows.size(), m.Ncols());
        int masked_row = 1;
        for (int r = 1; r <= m.Nrows(); r++)
        {
            if (std::find(masked_rows.begin(), masked_rows.end(), r) == masked_rows.end())
            {
                masked.Row(masked_row) = m.Row(r);
                masked_row++;
            }
        }
        return masked;
    }
}

ReturnMatrix MaskRows(ColumnVector v, vector<int> masked_rows)
{
    if (masked_rows.size() == 0)
    {
        return v;
    }
    else
    {
        NEWMAT::ColumnVector masked(v.Nrows() - masked_rows.size());
        int masked_row = 1;
        for (int r = 1; r <= v.Nrows(); r++)
        {
            if (std::find(masked_rows.begin(), masked_rows.end(), r) == masked_rows.end())
            {
                masked(masked_row) = v(r);
                masked_row++;
            }
        }
        return masked;
    }
}
}

double gammaln(double x)
{
    ColumnVector series(7);
    series << 2.5066282746310005 << 76.18009172947146 << -86.50532032941677 << 24.01409824083091
           << -1.231739572450155 << 0.1208650973866179e-2 << -0.5395239384953e-5;

    double total = 1.000000000190015;
    for (int i = 2; i <= series.Nrows(); i++)
        total += series(i) / (x + i - 1);

    return log(series(1) * total / x) + (x + 0.5) * log(x + 5.5) - x - 5.5;
}

double DescendingZeroFinder::FindZero() const
{
    double lower = searchMin;
    double upper = searchMax;
    double atLower, atUpper;

    double atSearchGuess = fcn(searchGuess);

    if (verbosity >= 2)
        LOG_ERR("IG: f(" << (searchGuess) << ") == " << atSearchGuess << endl);

    if (atSearchGuess < 0)
    {
        upper = searchGuess;
        atUpper = atSearchGuess;
        atLower = fcn(lower);
        if (verbosity >= 2)
            LOG_ERR("LG: f(" << (lower) << ") == " << atLower << endl);
        if (atLower <= 0)
            return lower; // hit the limit
    }
    else
    {
        lower = searchGuess;
        atLower = atSearchGuess;
        atUpper = fcn(upper);
        if (verbosity >= 2)
            LOG_ERR("UG: f(" << (upper) << ") == " << atUpper << endl);
        if (atUpper >= 0)
            return upper; // hit the limit
    }

    int evals = maxEvaluations - 2;
    double maxJump = searchScale;
    double guess, prevGuess = searchGuess;
    //    double nonlinearity = 10; // force a bisection first time

    // Interpolation only
    while ((evals > 0) && (upper - lower > tolX || atLower - atUpper > tolY
                              || upper / lower > ratioTolX || atUpper / atLower > ratioTolY))
    {
        guess = guesstimator->GetGuess(lower, upper, atLower, atUpper);

        if (lower == guess || guess == upper)
        {
            // This should only happen if we're near the limits of
            // double precision (constant factor probably depends on
            // the guesstimator)
            assert(upper - lower <= 2 * fabs(lower) * numeric_limits<double>::epsilon());
            LOG_ERR("DescendingZeroFinder: giving up without reaching tolerances because we're at "
                    "the limits of double precision!\n");
            LOG_ERR("Lower: f(" << lower << ") == " << atLower << endl
                                << "Upper: f(" << upper << ") == " << atUpper << endl);
            break;
        }

        assert(lower < guess && guess < upper);
        // cout << lower << "<" << guess << "<" << upper << endl;
        if (!fcn.PickFasterGuess(&guess, lower, upper))
            evals--; // only count the non-cached evaluations.
        // cout << lower << "<" << guess << "<" << upper << endl;
        assert(lower < guess && guess < upper);

        if (guess - prevGuess > maxJump)
            guess = prevGuess + maxJump;
        else if (guess - prevGuess < -maxJump)
            guess = prevGuess - maxJump;

        maxJump *= searchScaleGrowth;
        prevGuess = guess;

        // never mind -- already out of bounds
        // from checking the limits above.

        //        cout << "upper" << "\t" << "guess" << "\t" << "lower" << endl;
        //        cout << upper << "\t" << guess << "\t" << lower << endl;
        assert(guess > lower);
        assert(guess < upper);
        double atGuess = fcn(guess);
        //        cout << atUpper << "\t" << atGuess << "\t" << atLower << endl << endl;

        if (verbosity >= 2)
            LOG_ERR("NG: f(" << (guess) << ") == " << atGuess << endl);

        //        double atGuessIfLinear =
        //            ( atLower+(atUpper-atLower)*(guess-lower)/(upper-lower) - atGuess );
        //        if (atGuess > atGuessIfLinear)
        //            nonlinearity = (atGuess-atGuessIfLinear)/(atGuess-atUpper);
        //        else
        //            nonlinearity = (atGuessIfLinear-atGuess)/(atLower-atGuess);

        //        cout << "Nonlinearity = " << nonlinearity << endl;
        //        cout << "atGuess = " << atGuess
        //             << "\natGuessIfLinear = " << atGuessIfLinear
        //             << "\natLower = " << atLower
        //             << "\natUpper = " << atUpper << endl;

        if (atGuess < 0)
        {
            upper = guess;
            atUpper = atGuess;
        }
        else
        {
            lower = guess;
            atLower = atGuess;
        }
    }

    /*
     // One final interpolation -- not necessary, we could pick anything
     // between lower and upper really.
     guess = guesstimator->GetGuess(lower, upper, atLower, atUpper);
     assert( lower <= guess && guess <= upper );
     fcn.PickFasterGuess(&guess, lower, upper, true);
     assert( lower <= guess && guess <= upper );
     */

    // Pick either lower or upper bound, depending on which is closer to zero
    assert(atLower >= 0 && -atUpper >= 0);
    if (atLower < -atUpper)
        guess = lower;
    else
        guess = upper;

    if (verbosity >= 1)
        LOG_ERR("Final upper/lower ratio: " << (upper / lower) << endl);

    return guess;
}

double RiddlersGuesstimator::GetGuess(double lower, double upper, double atLower, double atUpper)
{
    // equations below: from NRIC, section 9.2.  Simpler than Brent, slightly less reliable.

    if (halfDone)
    {
        // Phase two: fancy estimation.
        halfDone = false;
        double x3, fx3;

        if (x1 == lower)
        {
            x3 = upper;
            fx3 = atUpper;
        }
        else
        {
            assert(x2 == upper);
            x3 = lower;
            fx3 = atLower;
        }

        assert(x1 < x3 && x3 < x2);
        assert(fx2 < fx1);

        WARN_ONCE("Riddler's Method; No special cases!");

        if (false) // x3 != (x1 + x2)/2)
        {
            LOG_ERR("x3 == " << x3 << ", x1 == " << x1 << ", x2 = " << x2 << endl);
            LOG_ERR("x3 - (x1+x2)/2 == " << x3 - (x1 + x2) / 2 << endl);
            WARN_ALWAYS("Riddler's Method: x3 != (x1+x2)/2");
        }
        else if (true) //(fx2 < fx3 && fx3 < fx1)
        {
            double s = (fx1 - fx2 > 0) ? +1.0 : -1.0; // s = sign(fx1-fx2)
            double x4 = x3 + (x3 - x1) * s * fx3 / sqrt(fx3 * fx3 - fx1 * fx2);

            WARN_ALWAYS("Riddler's Method: phase two");

            assert(lower < x4 && x4 < upper);

            return x4;
        }
        else
        {
            WARN_ALWAYS("Riddler's Method cheat: dropping back to the bisection method!");
        }
    }

    halfDone = true;
    // Phase one: just pick the midpoint, but save the values for phase two
    x1 = lower;
    fx1 = atLower;
    x2 = upper;
    fx2 = atUpper;

    WARN_ALWAYS("Riddler's Method: phase one");

    double x3 = (x1 + x2) / 2;
    return x3;
}
