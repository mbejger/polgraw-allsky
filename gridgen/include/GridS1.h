///////////////////////////////////////
// Name:        GridS1.h
// Author:      Andrzej Pisarski
// Copyright:   Andrzej Pisarski
// License:     CC-BY-NC-ND
// Created:     10/06/2015
// Modification:07/09/2017 A.P.
///////////////////////////////////////

#ifndef GridS1_H
#define GridS1_H

#include "DensityS1.h"
#include "FisherRM.h"
#include <vector>

class GridS1 : public DensityS1
{
    public:
        GridS1(const FisherRM *const fm =nullptr, unsigned int spindown_=1);

        std::vector<double> grid(double c0=0.75, double xi=0., unsigned int nfft=1048576) const;
        std::vector<double> convert(double c0, double xi, const std::vector<double>& grid_prim);
        // convert vector from hyper-sphere space to hyper-ellipsoid space

        double density(double c0=0.75,  unsigned int nfft=1048576) const;
        //double density(double c0,  unsigned int nfft, unsigned int data_length) const; // this same defintion in DensityS1
        // Need ephemeris to set data length

    protected:
        const FisherRM *const m_fm;
};

#endif // GridS1_H
