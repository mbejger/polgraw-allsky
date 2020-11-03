///////////////////////////////////////
// Name:        GridS1DS.h
// Author:      Andrzej Pisarski
// Copyright:   Andrzej Pisarski
// License:     CC-BY-NC-ND
// Created:     25/10/2019
// Modification:25/10/2019 A.P.
///////////////////////////////////////

#ifndef GridS1DS_H
#define GridS1DS_H

#include "DensityS1DS.h"
#include "FisherRMDS.h"
#include <vector>

class GridS1DS : public DensityS1DS
{
    public:
        GridS1DS(const FisherRMDS *const fmds = nullptr, unsigned int spindown_=1);

        std::vector<double> grid(double c0=0.75, double xi=0., unsigned int nfft=1048576, unsigned int data_length=344656) const;
        std::vector<double> convert(double c0, double xi, const std::vector<double>& grid_prim);
        // convert vector from hyper-sphere space to hyper-ellipsoid space
        // nfft=1048576 == 2^20, data_length=344656 - two sidereal days, delta_t_sample = 0.5 s.

        //double density(double c0,  unsigned int nfft, unsigned int data_length) const;

    protected:
        const FisherRMDS *const m_fmds;
};

#endif // GridS1DS_H
