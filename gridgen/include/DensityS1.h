///////////////////////////////////////
// Name:        DensityS1.h
// Author:      Andrzej Pisarski
// Copyright:   Andrzej Pisarski
// License:     CC-BY-NC-ND
// Created:     10/06/2015
// Modification:12/01/2020 A.P.
///////////////////////////////////////

#ifndef DENSITYS1_H
#define DENSITYS1_H

#include <vector>

class DensityS1
{
    public:
        DensityS1(unsigned int spindown_=1);
        double density(const std::vector<double>& basis) const; // need n-th base vectors in  n-th D space with hyper-spheres
        double density(double c0=0.75,  unsigned int nfft=1048576, unsigned int data_length=344656) const; // Do not need ephemeris (to set data length)
        std::vector<double> grid_prim(double c0, unsigned int nfft, unsigned int data_length) const; // return basis in space
                                                                                    // with hyper-spheres
        /// 'density' means this same as 'thickness'
    protected:
        unsigned int m_spindown;
        unsigned int m_dim;                                 // An* dimension
        unsigned int sd2dim(unsigned int spindown=1) const;  // convert spindown to An* dimension (all sky searches)
        unsigned int dim2sd(unsigned int dimension=4) const; // convert An* dimension to spindown (all sky searches)
};

#endif // DENSITYS1_H
