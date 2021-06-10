/*  dist_mvn.h - MultiVariate Normal distribution class/structure

 Adrian Groves, FMRIB Image Analysis Group

 Copyright (C) 2007-2008 University of Oxford  */

/*  CCOPYRIGHT */

#pragma once

#include "easylog.h"
#include "rundata.h"

#include "armawrap/newmat.h"
#include <ostream>
#include <string>
#include <vector>

/**
 * Multivariate Normal Distribution
 *
 * For a given number of parameters, the distribution
 * is defined by:
 *    - The means, one for each parameter
 *    - The covariances, or equivalently their precisions,
 *      one for each distinct pair of parameters and
 *      irrespective of order
 */
class MVNDist : public Loggable
{
public:
    /**
     * Load a per-voxel vector of MVN distributions from existing voxel data
     *
     * @param mvns One MVN for each voxel
     * @param mvns MVN in the form of voxel data as written by MVNDist::Save
     */
    static void Load(std::vector<MVNDist *> &mvns, NEWMAT::Matrix &voxel_data, EasyLog *log);

    /**
     * Load a per-voxel vector of MVN distributions from run data
     *
     * This may load from a NIFTI file or from data explicitly provided
     * using SetVoxelData.
     *
     * @param mvns Vector which should initially be empty. If non-empty, all entries must be NULL.
     * Will
     *             be resized to one MVN per voxel and the contents will be allocated using new and
     *             must be freed after use.
     * @param key Voxel data key in the run data. Equivalent to the option used to specify the
     *            filename using the command line tool
     * @param data Options and voxel data
     */
    static void Load(
        std::vector<MVNDist *> &mvns, const std::string &key, FabberRunData &data, EasyLog *log);

    /**
     * Save a per-voxel vector of MVN distributions.
     *
     * This may save to a NIFTI file if the run has been configured to save
     * files, or alternatively will save the matrix to the run data where
     * it can be retrieved using GetVoxelData.
     *
     * The matrix will be an N+1xN+1 symmetric matrix as described in /ref Load,
     * however only the lower triangle data is specified. When saving to a NIFTI
     * file NIFTI_INTENT_SYMMATRIX will be set. The resulting data is in the same
     * format expected by Load.
     *
     * @param mvns One MVN for each voxel
     */
    static void Save(const vector<MVNDist *> &mvns, const string &filename, FabberRunData &data);

    /**
     * Default constructor
     *
     * Size will be undefined -- will be fixed by first SetPrecisions/SetCovariance. Any
     * attempt to get data before this will throw an exception.
     */
    MVNDist(EasyLog *log = NULL);

    /**
     * Create distribution of known size
     */
    MVNDist(int dim, EasyLog *log = 0);

    /**
     * Copy constructor
     */
    MVNDist(const MVNDist &from);

    /**
     * Concatentate two MVNs
     */
    MVNDist(const MVNDist &from1, const MVNDist &from2);

    /**
     * Create from matrix file (VEST or ASCII)
     */
    MVNDist(const string filename, EasyLog *log = 0);

    /**
     * Copy using a subset of another MVN distribution's parameters
     */
    void CopyFromSubmatrix(const MVNDist &from, int first, int last, bool checkIndependence = true);

    /**
     * Get a subset of this MVN distribution as another MVN distribution
     */
    MVNDist GetSubmatrix(int first, int last, bool checkIndependence = true);

    /**
     * Assignment operator
     */
    MVNDist &operator=(const MVNDist &from);

    /**
     * Set the size
     *
     * Will resize means and covariances/precisions to match
     */
    void SetSize(int dim);

    /**
     * Get the size (number of parameters)
     *
     * @return size or -1 if not initialized
     */
    int GetSize() const;

    /**
     * Mean values of each parameter
     *
     * FIXME this is a public member which shouldn't be
     * resized - should encapsulate it to protect
     * against misuse
     */
    NEWMAT::ColumnVector means;

    /**
     * Get the precisions
     *
     * This is a symmetric matrix, i.e. the precision of two
     * parameters is unaffected by what order you specify
     * them in.
     *
     * Precision is the inverse of covariance, just
     * as in 1-dimension the precision of a random
     * variable is the reciprocol of the variance
     *
     * The reference returned is only valid temporarily -
     * it could become out of date if a subsequent call
     * to SetXXX is made. It should not be stored for future use.
     */
    const NEWMAT::SymmetricMatrix &GetPrecisions() const;

    /**
     * Get the covariances
     *
     * @see GetPrecisions
     */
    const NEWMAT::SymmetricMatrix &GetCovariance() const;

    /**
     * Set the precisions
     *
     * Also checks that size of matrix
     * is correct for current number of parameters
     *
     * Covariances will be updated lazily on next
     * call to GetCovariances
     */
    void SetPrecisions(const NEWMAT::SymmetricMatrix &from);

    /**
     * Set the covariances
     *
     * Also checks that size of matrix
     * is correct for current number of parameters
     *
     * Precisions will be updated lazily on next
     * call to GetPrecisions
     */
    void SetCovariance(const NEWMAT::SymmetricMatrix &from);

    /**
     * Load from matrix file
     *
     * The matrix may be in VEST or ASCII format
     *
     * The file must contain an N+1 x N+1 symmetric matrix.
     * This contains the covariances and the means in the
     * following example form (C=covariance, M=mean)
     *
     * CCCM
     * CCCM
     * CCCM
     * MMM1
     *
     * Precisions will be updated lazily on next
     * call to GetPrecisions
     *
     * @param filename VEST or ASCII matrix file
     */
    void LoadFromMatrix(const std::string &filename);

    /**
     * Dump info to output ostream
     */
    void Dump(std::ostream &os) const;

private:
    int m_size; // should only be changed explicitly

    // Mutable, because they're a lazily calculated on first use
    // even for const instance.
    mutable NEWMAT::SymmetricMatrix precisions;
    mutable NEWMAT::SymmetricMatrix covariance;
    mutable bool precisionsValid;
    mutable bool covarianceValid;
};

inline std::ostream &operator<<(std::ostream &out, const MVNDist &dist)
{
    dist.Dump(out);
    return out;
}
