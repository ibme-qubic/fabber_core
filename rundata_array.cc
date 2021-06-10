/**
 * fabber_rundata_array.cc
 *
 * Martin Craig
 *
 * Copyright (C) 2016 University of Oxford  */

/*  CCOPYRIGHT */

#include "rundata_array.h"

#include "easylog.h"
#include "rundata.h"

#include "armawrap/newmat.h"

#include <string>
#include <vector>

using namespace std;
using namespace NEWMAT;

void FabberRunDataArray::SetExtent(int nx, int ny, int nz, const int *mask)
{
    assert(nx > 0);
    assert(ny > 0);
    assert(nz > 0);
    assert(mask);

    FabberRunData::SetExtent(nx, ny, nz);

    int nv = nx * ny * nz;
    if (mask)
    {
        m_mask.insert(m_mask.end(), mask, mask + nv);
    }
    else
    {
        m_mask.resize(nv, 1);
    }

    int *maskPtr = &m_mask[0];
    int v = 0;
    Matrix coords(3, nv);
    for (int z = 0; z < nz; z++)
    {
        for (int y = 0; y < ny; y++)
        {
            for (int x = 0; x < nx; x++)
            {
                bool masked = (*maskPtr == 0);
                if (!masked)
                {
                    coords(1, v + 1) = x;
                    coords(2, v + 1) = y;
                    coords(3, v + 1) = z;
                    ++v;
                }
                ++maskPtr;
            }
        }
    }

    coords = coords.Columns(1, v);
    SetVoxelCoords(coords);
}

void FabberRunDataArray::GetVoxelDataArray(string key, float *data)
{
    assert(data);
    Matrix mdata = FabberRunData::GetVoxelData(key);
    int nt = mdata.Nrows();
    float *dataPtr = data;

    for (int t = 0; t < nt; t++)
    {
        int v = 0;
        const int *maskPtr = &m_mask[0];
        for (int z = 0; z < m_extent[2]; z++)
        {
            for (int y = 0; y < m_extent[1]; y++)
            {
                for (int x = 0; x < m_extent[0]; x++)
                {
                    bool masked = (*maskPtr == 0);
                    if (!masked)
                        *dataPtr = mdata(t + 1, v + 1);
                    else
                        *dataPtr = 0;
                    ++dataPtr;
                    ++maskPtr;
                    if (!masked)
                        ++v;
                }
            }
        }
    }
}

void FabberRunDataArray::SetVoxelDataArray(string key, int data_size, const float *data)
{
    assert(data);
    const float *dataPtr = data;
    int num_voxels = m_extent[0] * m_extent[1] * m_extent[2];
    Matrix matrixData(data_size, num_voxels);

    int v = 0;
    for (int t = 0; t < data_size; t++)
    {
        v = 0;
        const int *maskPtr = &m_mask[0];
        for (int z = 0; z < m_extent[2]; z++)
        {
            for (int y = 0; y < m_extent[1]; y++)
            {
                for (int x = 0; x < m_extent[0]; x++)
                {
                    bool masked = (*maskPtr == 0);
                    if (!masked)
                    {
                        matrixData(t + 1, v + 1) = *dataPtr;
                    }
                    ++dataPtr;
                    ++maskPtr;
                    if (!masked)
                        ++v;
                }
            }
        }
    }
    matrixData = matrixData.Columns(1, v);
    FabberRunData::SetVoxelData(key, matrixData);
}
