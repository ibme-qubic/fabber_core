/**
 * fabber_rundata_newimage.h
 *
 * Extends FabberRunData with loading and saving of data to
 * NIFTI files. This is the only part of Fabber which depends
 * on NEWIMAGE and can be disabled if compiling for an environment
 * which has its own image load/save facilities.
 *
 * Martin Craig
 *
 * Copyright (C) 2016 University of Oxford  */

/*  CCOPYRIGHT */

#pragma once

#ifndef NO_NEWIMAGE
#include "rundata.h"

#include "newimage/newimage.h"
#include "armawrap/newmat.h"

#include <string>

/**
 * Run data which uses NEWIMAGE to load NIFTII files
 */
class FabberRunDataNewimage : public FabberRunData
{
public:
    FabberRunDataNewimage(bool compat_options = true);

    void SetExtentFromData();
    const NEWMAT::Matrix &LoadVoxelData(const std::string &filename);
    virtual void SaveVoxelData(
        const std::string &filename, NEWMAT::Matrix &data, VoxelDataType data_type = VDT_SCALAR);

private:
    void SetCoordsFromExtent(int nx, int ny, int nz);
    NEWIMAGE::volume<float> m_mask;
    bool m_have_mask;
};

#endif /* NO_NEWIMAGE */
