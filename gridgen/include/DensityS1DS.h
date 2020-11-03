///////////////////////////////////////
// Name:        DensityS1DS.h
// Author:      Andrzej Pisarski
// Copyright:   Andrzej Pisarski
// License:     CC-BY-NC-ND
// Created:     30/07/2019
// Modification:30/07/2019 A.P.
///////////////////////////////////////

#ifndef DENSITYS1DS_H
#define DENSITYS1DS_H

#include <vector>

class DensityS1DS
{
    public:
        DensityS1DS(unsigned int spindown_=1);
        double density(const std::vector<double>& basis) const; // need n-th base vectors in  n-th D space with hyper-spheres
        double density(double c0=0.75,  unsigned int nfft=1048576, unsigned int data_length=344656) const; // Do not need ephemeris (to set data length)
        std::vector<double> grid_prim(double c0, unsigned int nfft, unsigned int data_length) const; // return basis in space
                                                                                    // with hyper-spheres
    protected:
        unsigned int m_spindown;
        unsigned int m_dim;                                    // An* dimension
        unsigned int sd2dim(unsigned int spindown=1) const;  // convert spindown to An* dimension (directed searches)
        unsigned int dim2sd(unsigned int dimension=2) const; // convert An* dimension to spindown (directed searches)
};

#endif // DENSITYS1DS_H
