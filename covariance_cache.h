#pragma once
/*  covariance_cache.h - Part of the spatial VB framework

 Adrian Groves & Michael Chappell, FMRIB Image Analysis Group & IBME QuBIc Group

 Copyright (C) 2007-2015 University of Oxford  */

/*  CCOPYRIGHT */

#include "easylog.h"

#include "armawrap/newmat.h"

#include <map>
#include <string>
#include <utility>

class CovarianceCache : public Loggable
{
public:
    void CalcDistances(const NEWMAT::Matrix &voxelCoords, const std::string &distanceMeasure);
    const NEWMAT::SymmetricMatrix &GetDistances() const;
    const NEWMAT::ReturnMatrix GetC(double delta) const; // quick to calculate
    const NEWMAT::SymmetricMatrix &GetCinv(double delta) const;
    const NEWMAT::SymmetricMatrix &GetCiCodistCi(double delta, double *CiCodistTrace = NULL) const;

    /**
     * If there's a cached value in (lower, upper), set *guess = value and
     * return true; otherwise return false and don't change *guess.
     */
    bool GetCachedInRange(
        double *guess, double lower, double upper, bool allowEndpoints = false) const;

private:
    typedef std::map<double, NEWMAT::SymmetricMatrix> Cinv_cache_type;
    typedef std::map<double, std::pair<NEWMAT::SymmetricMatrix, double> > CiCodistCi_cache_type;

    NEWMAT::SymmetricMatrix m_distances;
    mutable Cinv_cache_type m_cinv_cache;
    mutable NEWMAT::SymmetricMatrix m_cinv;
    mutable CiCodistCi_cache_type CiCodistCi_cache;
};
