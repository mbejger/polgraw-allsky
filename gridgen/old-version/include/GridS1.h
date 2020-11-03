///////////////////////////////////////
// Name:        GridS1.h
// Author:      Andrzej Pisarski
// Copyright:   Andrzej Pisarski
// License:     CC-BY-NC-ND
// Created:     10/06/2015
///////////////////////////////////////

#ifndef GridS1_H
#define GridS1_H

#include "DensityS1.h"
#include "FisherRM.h"
#include <vector>

class GridS1 : public DensityS1
{
    public:
        GridS1(const FisherRM *const);

        std::vector<double> grid(double, double, unsigned int) const;
        std::vector<double> convert(double, double, const std::vector<double>&); // convert vector
                                                    // from hyper-sphere space to hyper-ellipsoid space
        double density(double, unsigned int) const;
        double density(double, unsigned int, unsigned int) const;     // Need ephemeris to set data length

    protected:
        const FisherRM *const m_fm;
};

#endif // GridS1_H
