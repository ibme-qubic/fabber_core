#include "covariance_cache.h"

#include "easylog.h"
#include "rundata.h"

#include "armawrap/newmat.h"

#include <map>
#include <math.h>
#include <string>
#include <utility>

using namespace std;
using namespace NEWMAT;

// Caching is disabled currently
#define NOCACHE 1

// Euclidian distance
static double dist_euclid(double dx, double dy, double dz)
{
    return sqrt(dx * dx + dy * dy + dz * dz);
}

// Manhattan distance
static double dist_manh(double dx, double dy, double dz)
{
    return fabs(dx) + fabs(dy) + fabs(dz);
}
// Almost-squared Euclidian distance
static double dist_sqeuclid(double dx, double dy, double dz)
{
    return pow(dx * dx + dy * dy + dz * dz, 0.995);
}

/**
 * Calculate a distance matrix
 *
 * FIXME voxelCoords should really be in MM, not indices; only really matters
 * if it's aniostropic or you're using the smoothness values directly.
 *
 * @param voxelCoords List of voxel co-ordinates as a matrix: column = voxel
 * @param distanceMeasure How to measure distance:
 *   dist1 = Euclidian distance,
 *   dist2 = squared Euclidian distance,
 *   mdist = Manhattan distance (|dx| + |dy|)
 */
void CovarianceCache::CalcDistances(
    const NEWMAT::Matrix &voxelCoords, const string &distanceMeasure)
{
    const int nVoxels = voxelCoords.Ncols();

    // Create 3 column vectors, one to hold X co-ordinates, one Y and one Z
    ColumnVector positions[3];
    positions[0] = voxelCoords.Row(1).t();
    positions[1] = voxelCoords.Row(2).t();
    positions[2] = voxelCoords.Row(3).t();

    // dimSize is already included in voxelCoords
    // FIXME not obvious that it is, if not this should
    // be the dimensions of a voxel in mm
    const double dimSize[3] = { 1.0, 1.0, 1.0 };

    if (nVoxels > 7500)
    {
        LOG << "SpatialVariationalBayes::CalcDistances Over "
            << int(2.5 * nVoxels * nVoxels * 8 / 1e9)
            << " GB of memory will be used just to calculate "
            << "the distance matrix.  Hope you're not trying to invert this "
               "sucker!\n"
            << endl;
    }

    // 3 NxN symmetric matrices where N is the number of voxels
    // Each entry gives the absolute difference between
    // X/Y/Z co-ordinates respectively (FIXME in millimetres
    // supposedly but perhaps not given comments above)
    SymmetricMatrix relativePos[3];

    // Column vector with one entry per voxel, all set to 1 initially
    ColumnVector allOnes(nVoxels);
    allOnes = 1.0;

    // FIXME do not understand why this works!
    for (int dim = 0; dim < 3; dim++)
    {
        Matrix rel = dimSize[dim] * (positions[dim] * allOnes.t() - allOnes * positions[dim].t());
        assert(rel == -rel.t());
        assert(rel.Nrows() == nVoxels);
        // Down-convert to symmetric matrix (lower triangle so all positive?)
        relativePos[dim] << rel;
    }

    // Distances is an NxN symmetric matrix where N is number of voxels
    m_distances.ReSize(nVoxels);

    // Clunky but only run once so performance probably not significant
    for (int a = 1; a <= nVoxels; a++)
    {
        for (int b = 1; b <= a; b++)
        {
            if (distanceMeasure == "dist1")
                m_distances(a, b)
                    = dist_euclid(relativePos[0](a, b), relativePos[1](a, b), relativePos[2](a, b));
            else if (distanceMeasure == "dist2")
                m_distances(a, b) = dist_sqeuclid(
                    relativePos[0](a, b), relativePos[1](a, b), relativePos[2](a, b));
            else if (distanceMeasure == "mdist")
                m_distances(a, b)
                    = dist_manh(relativePos[0](a, b), relativePos[1](a, b), relativePos[2](a, b));
            else
                throw InvalidOptionValue(
                    "distance-measure", distanceMeasure, "Unrecognized distance measure");
        }
    }
}

const NEWMAT::SymmetricMatrix &CovarianceCache::GetDistances() const
{
    return m_distances;
}
const ReturnMatrix CovarianceCache::GetC(double delta) const
{
    const int Nvoxels = m_distances.Nrows();

    SymmetricMatrix C(Nvoxels);
    if (delta == 0)
    {
        C = IdentityMatrix(Nvoxels);
    }
    else
    {
        for (int a = 1; a <= Nvoxels; a++)
            for (int b = 1; b <= a; b++)
                C(a, b) = exp(-0.5 * m_distances(a, b) / delta);
    }

    // NOTE: when m_distances = squared distance, prior is equivalent to white
    // noise smoothed with a Gaussian with sigma^2 = 2*delta (haven't actually
    // double-checked this yet).
    // BEWARE: delta is measured in millimeters!! (based on NIFTI file info).

    C.Release();
    return C;
}

bool CovarianceCache::GetCachedInRange(
    double *guess, double lower, double upper, bool allowEndpoints) const
{
    assert(guess != NULL);
    const double initialGuess = *guess;
    if (!(lower < initialGuess && initialGuess < upper))
    {
        LOG << "SpatialVariationalBayes::Uh-oh... lower = " << lower
            << ", initialGuess = " << initialGuess << ", upper = " << upper << endl;
    }
    assert(lower < initialGuess && initialGuess < upper);

    Cinv_cache_type::iterator it = m_cinv_cache.lower_bound(lower);
    if (it == m_cinv_cache.end())
        return false;
    if (it->first == lower && !allowEndpoints)
        ++it;
    if (it == m_cinv_cache.end())
        return false;
    if (it->first > upper)
        return false;
    if (it->first == upper && !allowEndpoints)
        return false;

    // Success -- we have at least one fast guess!
    *guess = it->first;

    // Can we find a better one?
    while (++it != m_cinv_cache.end() && it->first <= upper)
    {
        if (it->first == upper && !allowEndpoints)
            break;

        if (it->first < initialGuess || it->first - initialGuess < initialGuess - *guess)
            *guess = it->first;
    }

    assert(lower < *guess && *guess < upper);

    return true;
}

const SymmetricMatrix &CovarianceCache::GetCinv(double delta) const
{
#ifdef NOCACHE
    WARN_ONCE("CovarianceCache::GetCinv Cache is disabled to avoid memory problems!");
    m_cinv = GetC(delta);

    // Appears to be memory bug when inverting a 0x0 matrix
    if (m_cinv.Nrows() > 0)
    {
        m_cinv = m_cinv.i();
    }

    return m_cinv;
#else
    if (m_cinv_cache[delta].Nrows() == 0)
    {
        m_cinv_cache[delta] = GetC(delta).i();
    }

    return m_cinv_cache[delta];
#endif
}

const SymmetricMatrix &CovarianceCache::GetCiCodistCi(double delta, double *CiCodistTrace) const
{
    if (CiCodistCi_cache[delta].first.Nrows() == 0)
    {
#ifdef NOCACHE
        CiCodistCi_cache.clear();
#endif
        GetCinv(delta); // for sensible messages, make sure cache hits

        Matrix CiCodist = GetCinv(delta) * SP(GetC(delta), m_distances);
        CiCodistCi_cache[delta].second = CiCodist.Trace();
        Matrix CiCodistCi_tmp = CiCodist * GetCinv(delta);
        CiCodistCi_cache[delta].first << CiCodistCi_tmp; // Force symmetric

        { // check something
            double maxAbsErr
                = (CiCodistCi_cache[delta].first - CiCodistCi_tmp).MaximumAbsoluteValue();
            if (maxAbsErr > CiCodistCi_tmp.MaximumAbsoluteValue() * 1e-5)
            // If that test fails, you're probably in trouble.
            // Reducing it to e.g. 1e-5 (to make dist2 work)
            //   => non-finite alpha in iteration 2
            // Reduced it to 1e-5 to get mdist to work....
            //   => same result.  (oops, was mdist2)
            // Reduced it to 1e-5 to get mdist to work, again...
            //   =>
            {
                LOG << "CovarianceCache::GetCiCodistCi matrix not symmetric!\nError = " << maxAbsErr
                    << ", maxabsvalue = " << CiCodistCi_tmp.MaximumAbsoluteValue() << endl;
                assert(false);
            }
        }
    }

    if (CiCodistTrace != NULL)
        (*CiCodistTrace) = CiCodistCi_cache[delta].second;
    return CiCodistCi_cache[delta].first;
}
