/*  dist_mvn.cc - MultiVariate Normal distribution class/structure

 Adrian Groves, FMRIB Image Analysis Group

 Copyright (C) 2007-2008 University of Oxford  */

/*  CCOPYRIGHT */

#include "dist_mvn.h"
#include "easylog.h"
#include "tools.h"

#include <math.h>
#include "armawrap/newmat.h"

using fabber::read_matrix_file;
using namespace NEWMAT;

// Constructors

MVNDist::MVNDist(EasyLog *log)
    : Loggable(log)
    , m_size(-1)
    , precisionsValid(false)
    , covarianceValid(false)
{
}

MVNDist::MVNDist(int dim, EasyLog *log)
    : Loggable(log)
    , m_size(-1)
    , precisionsValid(false)
    , covarianceValid(false)
{
    SetSize(dim);
}

MVNDist::MVNDist(const MVNDist &from)
    : Loggable(from.m_log)
    , m_size(-1)
    , precisionsValid(false)
    , covarianceValid(false)
{
    *this = from;
}

MVNDist::MVNDist(const string filename, EasyLog *log)
    : Loggable(log)
    , m_size(-1)
    , precisionsValid(false)
    , covarianceValid(false)
{
    LoadFromMatrix(filename);
}

MVNDist::MVNDist(const MVNDist &from1, const MVNDist &from2)
    : Loggable(from1.m_log)
    , m_size(-1)
    , precisionsValid(false)
    , covarianceValid(false)
{
    SetSize(from1.m_size + from2.m_size);

    means = from1.means & from2.means;

    // Always duplicate the covariances (even if this means some recalculation)
    // Otherwise if we use precisions.i(), zeros won't stay exactly zero
    covariance = 0;
    for (int row=1; row <= from1.m_size; row++)
    {
        for (int col=1; col <= from1.m_size; col++)
        {
            try
            {
                covariance(row, col) = from1.GetCovariance()(row, col);
            }
            catch (Exception)
            {
                covariance(row, col) = 0;
            }
        }
    }
    for (int row=1; row <= from2.m_size; row++)
    {
        for (int col=1; col <= from2.m_size; col++)
        {
            try
            {
                covariance(from1.m_size + row, from1.m_size + col) = from2.GetCovariance()(row, col);
            }
            catch (Exception)
            {
                covariance(row, col) = 0;
            }
        }
    }

    assert(means.Nrows() == m_size);
}

MVNDist &MVNDist::operator=(const MVNDist &from)
{
    // Special case: assignment to self (is a no-op)
    if (&from == this)
        return *this;

    m_log = from.m_log;

    // Special case: assigned from an uninitialized MVNDist
    if (from.m_size == -1)
    {
        m_size = -1;
        precisionsValid = covarianceValid = false;
        // Free any memory from previous instance
        SetSize(1);
        return *this;
    }

    assert(from.m_size == from.means.Nrows());

    SetSize(from.m_size);
    means = from.means;
    precisionsValid = from.precisionsValid;
    covarianceValid = from.covarianceValid;
    if (precisionsValid)
        precisions = from.precisions;

    if (covarianceValid)
        covariance = from.covariance;

    assert(means.Nrows() == m_size);
    return *this;
}

MVNDist MVNDist::GetSubmatrix(int first, int last, bool checkIndependence)
{
    MVNDist ret;
    ret.CopyFromSubmatrix(*this, first, last, checkIndependence);
    return ret;
}

void MVNDist::CopyFromSubmatrix(const MVNDist &from, int first, int last, bool checkIndependence)
{
    SetSize(last - first + 1);
    means = from.means.Rows(first, last);
    precisionsValid = from.precisionsValid;
    covarianceValid = from.covarianceValid;
    if (precisionsValid)
        precisions = from.precisions.SymSubMatrix(first, last);

    if (covarianceValid)
        covariance = from.covariance.SymSubMatrix(first, last);

    assert(means.Nrows() == m_size);

    if (checkIndependence)
    {
        Matrix deps1 = from.GetCovariance().Rows(first, last).Columns(1, first - 1);
        Matrix deps2
            = from.GetCovariance().Rows(first, last).Columns(last + 1, from.covariance.Ncols());
        if (!deps1.IsZero() || !deps2.IsZero())
            throw FabberRunDataError(
                "Covariance found in part of MVN that should be independent from the rest!");
    }
}

int MVNDist::GetSize() const
{
    assert(m_size == means.Nrows() || m_size < 0);
    return m_size;
}

void MVNDist::SetSize(int dim)
{
    if (dim <= 0)
        throw FabberInternalError("MVNDist::SetSize dim<=0");

    assert(means.Nrows() == m_size || m_size < 0);

    if (m_size != dim)
    {
        m_size = dim;
        means.ReSize(dim);
        means = 0;
        precisions = IdentityMatrix(dim);
        covariance = IdentityMatrix(dim);
    }
    precisionsValid = true;
    covarianceValid = true;

    assert(means.Nrows() == m_size);
    assert(precisions.Nrows() == m_size);
    assert(covariance.Nrows() == m_size);
}

const SymmetricMatrix &MVNDist::GetPrecisions() const
{
    if (m_size == -1)
        throw FabberInternalError("MVNDist::GetPrecisions size = -1 (uninitialized)");
    assert(means.Nrows() == m_size);

    // If the covariances have been changed then the
    // precisions are out of date and need to be
    // recalculated
    if (!precisionsValid)
    {
        assert(covarianceValid);
        // precisions and precisionsValid are mutable,
        // so we can change them even in a const function
        try
        {
            precisions = covariance.i();
        }
        catch (Exception)
        {
            // Failure to invert matrix - this hack adds a tiny amount to the diagonal and tries
            // again
            WARN_ONCE("MVN precision (m_size==" + stringify(m_size)
                + ") was singular, adding 1e-10 to diagonal");
            LOG << means.t() << endl;
            LOG << covariance << endl;
            precisions = (covariance + IdentityMatrix(m_size) * 1e-10).i();
        }
        precisionsValid = true;
    }
    assert(means.Nrows() == m_size);
    assert(precisions.Nrows() == m_size);
    return precisions;
}

const SymmetricMatrix &MVNDist::GetCovariance() const
{
    if (m_size == -1)
        throw FabberInternalError("MVNDist::GetCovariance size = -1 (uninitialized)");
    assert(means.Nrows() == m_size);

    // If the precisions have been changed then the
    // covariances are out of date and need to be
    // recalculated
    if (!covarianceValid)
    {
        assert(precisionsValid);
        // covariance and covarianceValid are mutable,
        // so we can change them even in a const function
        try
        {
            covariance = precisions.i();
        }
        catch (Exception)
        {
            // Failure to invert matrix - this hack adds a tiny amount to the diagonal and tries
            // again
            WARN_ONCE("MVN precision (m_size==" + stringify(m_size)
                + ") was singular, adding 1e-10 to diagonal");
            LOG << means.t() << endl;
            LOG << precisions << endl;
            covariance = (precisions + IdentityMatrix(m_size) * 1e-10).i();
        }
        covarianceValid = true;
    }
    assert(means.Nrows() == m_size);
    assert(covariance.Nrows() == m_size);
    return covariance;
}

void MVNDist::SetPrecisions(const SymmetricMatrix &from)
{
    assert(from.Nrows() == m_size);
    assert(means.Nrows() == m_size);
    precisions = from;
    precisionsValid = true;
    covarianceValid = false;
    assert(means.Nrows() == m_size);
}

void MVNDist::SetCovariance(const SymmetricMatrix &from)
{
    assert(from.Nrows() == m_size);
    assert(means.Nrows() == m_size);
    covariance = from;
    covarianceValid = true;
    precisionsValid = false;
    assert(means.Nrows() == m_size);
}

void MVNDist::LoadFromMatrix(const string &filename)
{
    LOG << "MVNDist::Reading MVN from file '" << filename << "'...\n";
    Matrix mat = read_matrix_file(filename);

    // Format: [covariance means(:); means(:)' 1.0]
    const int N = mat.Nrows() - 1;

    if (N < 1 || mat != mat.t() || mat(N + 1, N + 1) != 1.0)
    {
        LOG << "N == " << N << ", matrix:\n" << mat;
        throw InvalidOptionValue(filename, "",
            "MVNs must be symmetric matrices (format = [covariance means(:); means(:) 1.0])");
    }
    SetSize(N);
    means = mat.Column(m_size + 1).Rows(1, N);

    SymmetricMatrix sym;
    sym << mat.SubMatrix(1, N, 1, N);
    SetCovariance(sym);

    assert(means.Nrows() == m_size);
}

void MVNDist::Load(
    vector<MVNDist *> &mvns, const string &filename, FabberRunData &data, EasyLog *log)
{
    // Input matrix contains 3d voxels with the
    // 4th dimension containing the covariances
    // and means in the format as described in
    // Load. First this is converted into
    // a matrix whose columns are the voxels
    // and rows are the data
    Matrix voxel_data = data.GetVoxelData(filename);
    MVNDist::Load(mvns, voxel_data, log);
}

void MVNDist::Load(vector<MVNDist *> &mvns, Matrix &voxel_data, EasyLog *log)
{
    // Prepare an output vector of the correct size
    const int nVoxels = voxel_data.Ncols();
    if (nVoxels == 0)
    {
        throw FabberRunDataError("MVNDist::Load - Voxel data is empty");
    }

    for (unsigned i = 0; i < mvns.size(); i++)
        assert(mvns[i] == NULL); // should've deleted everything first.
    mvns.resize(nVoxels, NULL);

    // This formula for the number of parameters, P given
    // the number of rows, N, is found by inverting the formula
    // for N as a function of P given in Load, using the quadratic
    // formula.
    const int nParams = ((int)sqrt(double(8 * voxel_data.Nrows() + 1)) - 3) / 2;
    if (voxel_data.Nrows() != nParams * (nParams + 1) / 2 + nParams + 1)
    {
        throw FabberRunDataError("MVNDist::Load  - Incorrect number of rows for an MVN input");
    }

    SymmetricMatrix tmp(nParams);

    // Create a new MVN dist for each voxel,
    // and set the covariances and the means from
    // the data in the symmetric matrix
    for (int vox = 1; vox <= nVoxels; vox++)
    {
        MVNDist *mvn = new MVNDist(nParams, log);

        int index = 0;
        for (int r = 1; r <= nParams; r++)
            for (int c = 1; c <= r; c++)
                tmp(r, c) = voxel_data(++index, vox);

        assert(index == nParams * (nParams + 1) / 2);
        mvn->SetCovariance(tmp);
        mvn->means = voxel_data.Column(vox).Rows(
            nParams * (nParams + 1) / 2 + 1, nParams * (nParams + 1) / 2 + nParams);

        if (voxel_data(voxel_data.Nrows(), vox) != 1)
        {
            throw FabberRunDataError(
                "MVNDist::Load - Voxel data does not contain a valid MVN - last value != 1");
        }
        assert(mvn->means.Nrows() == mvn->m_size);
        assert(mvns.at(vox - 1) == NULL);
        mvns[vox - 1] = mvn;
    }
}

void MVNDist::Save(const vector<MVNDist *> &mvns, const string &filename, FabberRunData &data)
{
    // Save the MVNs in a NIFTI file as a single NIFTI_INTENT_SYMMATRIX
    // last row/col is the means (1 in the corner).
    // Note that I'm using the 4th dim and should really be using the 5th,
    // according to the specification -- but I don't think it really matters.
    Matrix vols;

    const int nVoxels = mvns.size();
    int nParams = 0; // In case we have no voxels
    if (nVoxels != 0)
    {
        assert(nVoxels > 0 && mvns.at(0) != NULL);
        nParams = mvns.at(0)->means.Nrows();
    }

    // This is what the matrix will look like (C = covariances, M=means)
    //
    // CC..CM
    //  C..CM
    //   C.CM
    //    CCM
    //     CM
    //      1
    //
    // This explains the formula below to calculate the number of data
    // elements required. The other triangle of the matrix is inferred
    // because of NIFTI_INTENT_SYMMATRIX
    int nCov = nParams * (nParams + 1) / 2;
    vols.ReSize(nCov + nParams + 1, nVoxels);
    ColumnVector aOne(1);
    aOne = 1.0;

    for (int vox = 1; vox <= nVoxels; vox++)
    {
        // Each column contains the values a voxel
        // Covariances first, but only the lower triangular numbers
        // Then means, and finally the last 1 as required by format above
        //
        // Note that AsColumn for a SymmetricMatrix uses row ordering on the
        // lower triangular part, returning (1,1) (2,1) (2,2) (3,1).. as
        // required by NIFTI_INTENT_SYMMATRIX
        ColumnVector cov(nCov);
        int idx=1;
        for (int row=1; row <=nParams; row++)
        {
            for (int col=1; col<=row; col++)
            {
                cov(idx++) = mvns.at(vox - 1)->GetCovariance()(row, col);
            }
        }

        vols.Column(vox) = cov & mvns.at(vox - 1)->means & aOne;
    }

    data.SaveVoxelData(filename, vols, VDT_MVN);
}

void MVNDist::Dump(ostream &out) const
{
    out << "MVNDist, with m_size == " << m_size << ", precisionsValid == " << precisionsValid
        << ", covarianceValid == " << covarianceValid << endl;
    out << "  Means: " << means.t();
    if (precisionsValid || covarianceValid)
    {
        out << "  Covariance matrix:" << endl;
        for (int i = 1; i <= m_size; i++)
            out << "  " << GetCovariance().Row(i);
    }
    else
        out << "  Covariance undefined." << endl;

    assert(means.Nrows() == m_size);
}
