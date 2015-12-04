///////////////////////////////////////
// Name:        DensityS1.h
// Author:      Andrzej Pisarski
// Copyright:   Andrzej Pisarski
// License:     CC-BY-NC-ND
// Created:     10/06/2015
///////////////////////////////////////

#ifndef DENSITYS1_H
#define DENSITYS1_H

#include <vector>

class DensityS1
{
    public:
        DensityS1();
        double density(const std::vector<double>& ) const; // need 4 base vectors in  4D space with hyper-spheres
        double density(double, unsigned int, unsigned int) const; // Do not need ephemeris (to set data length)
        std::vector<double> grid_prim(double, unsigned int, unsigned int) const;    // return basis in space
                                                                                    // with hyper-spheres
};

#endif // DENSITYS1_H
