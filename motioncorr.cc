/*  fabber_mc.h - Motion correction class

 Adrian Groves & Michael Chappell, FMRIB Image Analysis Group & IBME QuBIc Group

 Copyright (C) 2007-2015 University of Oxford  */

/*  CCOPYRIGHT */

#ifdef __FABBER_MOTION

#include "motioncorr.h"

#include <newimage/newimageall.h>
#include <warpfns/warpfns.h>

#include "armawrap/newmat.h"

#include <cmath>
#include <stdexcept>
#include <vector>

using namespace NEWIMAGE;
using namespace NEWMAT;

MCobj::MCobj(FabberRunData &allData, int dof)
{
    Matrix coords = allData.GetVoxelCoords();
    std::vector<int> size(3);
    size[0] = coords.Row(1).Maximum() + 1;
    size[1] = coords.Row(2).Maximum() + 1;
    size[2] = coords.Row(3).Maximum() + 1;
    mask = volume<float>(size[0], size[1], size[2]);
    mask = 0.0;

    // cout << size[0] << ", " << size[1] << ", " << size[2] << endl;
    // cout << coords.Ncols() << endl;
    for (int v = 1; v <= coords.Ncols(); v++)
    {
        mask(coords(1, v), coords(2, v), coords(3, v)) = 1.0;
        //	cout << coords(1, v) << ", " << coords(2, v) << ", " <<  coords(3, v) << endl;
    }
    Matrix datamat = allData.GetMainVoxelData();
    wholeimage.setmatrix(datamat, mask);

    // the following sets up an initial zero deformation field
    userdof = dof;
    num_iter = 10;
    modelpred = wholeimage;
    modelpred = 0.0f;
    if (userdof > 12)
    {
        defx = modelpred;
        defx.setROIlimits(0, 2);
        defx.activateROI();
        defx = defx.ROI();
        defy = defx;
        defz = defx;
        // Unnecessary initialisations?!?
        tmpx = defx;
        tmpy = defx;
        tmpz = defx;
    }
    else
    {
        mcf.setparams("verbose", false);
    }
    affmat = IdentityMatrix(4);
    finalimage = modelpred;
}

void MCobj::run_mc(const Matrix &modelpred_mat, Matrix &finalimage_mat)
{
    modelpred.setmatrix(modelpred_mat, mask);
    if (userdof > 12)
    {
        UpdateDeformation(
            wholeimage, modelpred, num_iter, defx, defy, defz, finalimage, tmpx, tmpy, tmpz);
        defx = tmpx;
        defy = tmpy;
        defz = tmpz;
    }
    else
    {
        // mcf.register_volumes(4D reference,4D input image,refweight,inweight,4D output image);
        affmat = mcf.register_volumes(modelpred, wholeimage, mask, mask, finalimage);

        // apply transforms to wholeimage to get finalimage (the above is a dummy)
        for (int n = 0; n < wholeimage.maxt(); n++)
        {
            affine_transform(wholeimage[n], finalimage[n], affmat.Rows(n * 4 + 1, n * 4 + 4));
        }
    }
    finalimage_mat = finalimage.matrix(mask);
}

/////////////////// Diffeomorphic code

// assuming everything is in mm
void diffeomorphic_new(
    const volume4D<float> &input_velocity, volume4D<float> &output_def, int steps)
{
    float coeff = 1.0f / (2 ^ steps);
    // cerr << "DN::COEFF = " << coeff << endl;
    volume4D<float> prewarp;
    prewarp = input_velocity * coeff;

    convertwarp_rel2abs(prewarp);
    for (int i = 0; i < steps; i++)
    {
        // print_volume_info(prewarp,"DN::prewarp");
        concat_warps(prewarp, prewarp, output_def);
        // print_volume_info(output_def,"DN::output_def");
        prewarp = output_def;
    }
    convertwarp_abs2rel(output_def);
}

////////////////// Other code

void calculate_gradients(volume4D<float> &gradient_imagex, volume4D<float> &gradient_imagey,
    volume4D<float> &gradient_imagez, const volume4D<float> &wholeimage)
{
    if (gradient_imagex.tsize() != wholeimage.tsize())
        gradient_imagex = wholeimage;
    if (gradient_imagey.tsize() != wholeimage.tsize())
        gradient_imagey = wholeimage;
    if (gradient_imagez.tsize() != wholeimage.tsize())
        gradient_imagez = wholeimage;

    volume4D<float> gradient_all;

    for (int t = 0; t < wholeimage.tsize(); t++)
    {
        gradient(wholeimage[t], gradient_all);
        gradient_imagex[t] = gradient_all[0];
        gradient_imagey[t] = gradient_all[1];
        gradient_imagez[t] = gradient_all[2];
    }
}

void UpdateDeformation(const volume4D<float> &wholeimage, const volume4D<float> &modelpred,
    int no_iter, const volume4D<float> &prevdefx, const volume4D<float> &prevdefy,
    const volume4D<float> &prevdefz, volume4D<float> &finalimage, volume4D<float> &defx,
    volume4D<float> &defy, volume4D<float> &defz)
{
    int steps = 4;

    const float lamda = 10;
    std::cout << "The update was compiled!" << lamda << no_iter << std::endl;

    print_volume_info(modelpred, "modelpred");
    print_volume_info(wholeimage, "wholeimage");
    print_volume_info(prevdefx, "prevdefx");
    print_volume_info(prevdefy, "prevdefy");
    print_volume_info(prevdefz, "prevdefz");

    volume4D<float> H11(wholeimage);
    volume4D<float> H12(wholeimage);
    volume4D<float> H13(wholeimage);
    volume4D<float> H21(wholeimage);
    volume4D<float> H22(wholeimage);
    volume4D<float> H23(wholeimage);
    volume4D<float> H31(wholeimage);
    volume4D<float> H32(wholeimage);
    volume4D<float> H33(wholeimage);

    volume4D<float> AH11(wholeimage);
    volume4D<float> AH12(wholeimage);
    volume4D<float> AH13(wholeimage);
    volume4D<float> AH21(wholeimage);
    volume4D<float> AH22(wholeimage);
    volume4D<float> AH23(wholeimage);
    volume4D<float> AH31(wholeimage);
    volume4D<float> AH32(wholeimage);
    volume4D<float> AH33(wholeimage);

    volume4D<float> Det(wholeimage);
    volume4D<float> wholeimage1(wholeimage);
    volume4D<float> diffimage(wholeimage);
    volume4D<float> input_velocity;
    volume4D<float> gradient_imagex, gradient_imagey, gradient_imagez, gradient_all;
    input_velocity.addvolume(wholeimage[0]);
    input_velocity.addvolume(wholeimage[0]);
    input_velocity.addvolume(wholeimage[0]);
    volume4D<float> output_def;
    output_def.addvolume(wholeimage[0]);
    output_def.addvolume(wholeimage[0]);
    output_def.addvolume(wholeimage[0]);

    defx = prevdefx;
    defy = prevdefy;
    defz = prevdefz;
    H11 = 1;
    H12 = 0;
    H13 = 0;
    H21 = 0;
    H22 = 1;
    H23 = 0;
    H31 = 0;
    H32 = 0;
    H33 = 1;

    int xnum = wholeimage.xsize();
    int ynum = wholeimage.ysize();
    int znum = wholeimage.zsize();
    int tnum = wholeimage.tsize();

    for (int t = 0; t < tnum; t++)
    {
        input_velocity[0] = prevdefx[t];
        input_velocity[1] = prevdefy[t];
        input_velocity[2] = prevdefz[t];

        diffeomorphic_new(input_velocity, output_def, steps);

        std::cout << "diffeomorphic done" << lamda << std::endl;

        // apply_warp(const volume<T> & invol, volume<T> &outvol, const volume4D<float>&  warpvol)
        apply_warp(wholeimage[t], wholeimage1[t], output_def);
    }

    diffimage = wholeimage1 - modelpred;
    print_volume_info(diffimage, "diffimage");

    double sum = diffimage.sumsquares();
    double new_similarity = sum / (xnum * ynum * znum * tnum);

    calculate_gradients(gradient_imagex, gradient_imagey, gradient_imagez, wholeimage1);

    gradient_imagex = -gradient_imagex * diffimage;
    gradient_imagey = -gradient_imagey * diffimage;
    gradient_imagez = -gradient_imagez * diffimage;

    print_volume_info(gradient_imagex, "gradient_imagex");
    // print_volume_info(gradient_imagey,"gradient_imagey");
    // print_volume_info(gradient_imagez,"gradient_imagez");

    double diff_similarity = 1;
    int count = 0;

    while ((diff_similarity > 0) && ((count++) < no_iter))
    // while (((count++)<no_iter))
    {
        double old_similarity = new_similarity;

        defx += H11 * gradient_imagex + H12 * gradient_imagey + H13 * gradient_imagez;
        defx = smooth(defx, 2);
        print_volume_info(defx, "defx");
        defy += H21 * gradient_imagex + H22 * gradient_imagey + H23 * gradient_imagez;
        defy = smooth(defy, 2);
        defz += H31 * gradient_imagex + H32 * gradient_imagey + H33 * gradient_imagez;
        defz = smooth(defz, 2);

        for (int t = 0; t < tnum; t++)
        {
            input_velocity[0] = defx[t];
            input_velocity[1] = defy[t];
            input_velocity[2] = defz[t];
            // print_volume_info(input_velocity,"dn::input_velocity");

            diffeomorphic_new(input_velocity, output_def, steps);

            // print_volume_info(output_def,"output_def");
            apply_warp(wholeimage[t], wholeimage1[t], output_def);

            // gradient(wholeimage1[t],gradient_all);

            // gradient_imagex[t]= gradient_all[0];
            // gradient_imagex[t]=
            // gradient_imagex[t]/Max(gradient_imagex[t].max()-gradient_imagex[t].min(),1e-6);
            // gradient_imagex[t]= gradient_imagex[t];

            // gradient_imagey[t]= gradient_all[1];
            // gradient_imagey[t]=
            // gradient_imagex[t]/Max(gradient_imagey[t].max()-gradient_imagey[t].min(),1e-6);
            // gradient_imagey[t]= gradient_imagey[t];

            // gradient_imagez[t]= gradient_all[2];
            // gradient_imagez[t]=
            // gradient_imagez[t]/Max(gradient_imagez[t].max()-gradient_imagez[t].min(),1e-6);
            // gradient_imagez[t]= gradient_imagez[t];
        }

        calculate_gradients(gradient_imagex, gradient_imagey, gradient_imagez, wholeimage1);

        diffimage = wholeimage1 - modelpred;
        sum = diffimage.sumsquares();
        new_similarity = sum / (xnum * ynum * znum * tnum);
        diff_similarity = old_similarity - new_similarity;

        std::cout << "new_similarity=" << new_similarity << std::endl;

        gradient_imagex = -gradient_imagex * diffimage;
        gradient_imagey = -gradient_imagey * diffimage;
        gradient_imagez = -gradient_imagez * diffimage;

        // print_volume_info(gradient_imagex,"gradient_imagex");
        // print_volume_info(gradient_imagey,"gradient_imagey");
        // print_volume_info(gradient_imagez,"gradient_imagez");

        Det = (lamda + gradient_imagex * gradient_imagex)
                * ((lamda + gradient_imagey * gradient_imagey)
                          * (lamda + gradient_imagez * gradient_imagez)
                      - gradient_imagey * gradient_imagez * gradient_imagey * gradient_imagez)
            - gradient_imagex * gradient_imagey
                * (gradient_imagex * gradient_imagey * (lamda + gradient_imagez * gradient_imagez)
                      - gradient_imagex * gradient_imagez * gradient_imagey * gradient_imagez)
            + gradient_imagex * gradient_imagez
                * (gradient_imagex * gradient_imagey * gradient_imagey * gradient_imagez
                      - (lamda + gradient_imagey * gradient_imagey) * gradient_imagex
                          * gradient_imagez);

        print_volume_info(Det, "Det");

        AH11 = (lamda + gradient_imagex * gradient_imagex);
        AH12 = gradient_imagex * gradient_imagey;
        AH13 = gradient_imagex * gradient_imagez;
        AH21 = AH12;
        AH22 = (lamda + gradient_imagey * gradient_imagey);
        AH23 = gradient_imagey * gradient_imagez;
        AH31 = AH13;
        AH32 = AH23;
        AH33 = (lamda + gradient_imagez * gradient_imagez);

        H11 += (AH22 * AH33 - AH23 * AH32) / Det;
        H12 += (AH32 * AH13 - AH33 * AH12) / Det;
        H13 += (AH23 * AH12 - AH22 * AH13) / Det;
        H21 += (AH31 * AH23 - AH33 * AH21) / Det;
        H22 += (AH33 * AH11 - AH31 * AH13) / Det;
        H23 += (AH21 * AH13 - AH23 * AH11) / Det;
        H31 += (AH32 * AH21 - AH31 * AH22) / Det;
        H32 += (AH31 * AH12 - AH32 * AH11) / Det;
        H33 += (AH22 * AH11 - AH21 * AH12) / Det;
        print_volume_info(H11, "H11");
    }
    finalimage = wholeimage1;

    print_volume_info(finalimage, "finalimage");

    std::cout << "The update has finished" << lamda << std::endl;
}

// Concerns:
// things to check: is the usage of 'smooth' fine? Can you have the same input and output?
// Should the declarations be in a separate header file

// Things to remember (for myself): that here, updates are being scaled, while the total deformation
// field is being smoothed. Change if necessary.
// All of this in voxels- convert to mm by simply multiplying by voxel sizes.

#endif //__FABBER_MOTION
