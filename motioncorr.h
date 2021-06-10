/*  fabber_mc.h - Motion correction class

 Adrian Groves & Michael Chappell, FMRIB Image Analysis Group & IBME QuBIc Group

 Copyright (C) 2007-2015 University of Oxford  */

/*  CCOPYRIGHT */
#pragma once

#ifdef __FABBER_MOTION

#include "rundata.h"

#include <mcflirt/rigidreglib.h>
#include <newimage/newimage.h>
#include <warpfns/warpfns.h>

#include "armawrap/newmat.h"

//   NB: for now the mask should cover the *entire* image as we zero everything
//       outside of the mask, which is not good for registration
//       In future we'd need allData to be able to provide the original image (or something to)

class MCobj
{
public:
    MCobj(FabberRunData &allData, int dof);
    void run_mc(const NEWMAT::Matrix &modelpred_mat, NEWMAT::Matrix &finalimage_mat);
    void set_num_iter(int nit)
    {
        num_iter = nit;
    }

private:
    int userdof;  // anything over 13 is full nonlinear
    int num_iter; // default 10
    NEWIMAGE::volume<float> mask;
    NEWMAT::Matrix affmat;
    mcflirt mcf;
    NEWIMAGE::volume4D<float> defx;
    NEWIMAGE::volume4D<float> defy;
    NEWIMAGE::volume4D<float> defz;
    // things below are kept for efficiency (?) in order to avoid repeated allocation/destruction
    NEWIMAGE::volume4D<float> tmpx;
    NEWIMAGE::volume4D<float> tmpy;
    NEWIMAGE::volume4D<float> tmpz;
    NEWIMAGE::volume4D<float> modelpred;
    NEWIMAGE::volume4D<float> finalimage;
    NEWIMAGE::volume4D<float> wholeimage;
};

void UpdateDeformation(const NEWIMAGE::volume4D<float> &wholeimage,
    const NEWIMAGE::volume4D<float> &modelpred, int no_iter,
    const NEWIMAGE::volume4D<float> &prevdefx, const NEWIMAGE::volume4D<float> &prevdefy,
    const NEWIMAGE::volume4D<float> &prevdefz, NEWIMAGE::volume4D<float> &finalimage,
    NEWIMAGE::volume4D<float> &defx, NEWIMAGE::volume4D<float> &defy,
    NEWIMAGE::volume4D<float> &defz);

#endif //__FABBER_MOTION
