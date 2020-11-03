///////////////////////////////////////
// Name:        DensityS2.h
// Author:      Andrzej Pisarski
// Copyright:   Andrzej Pisarski
// License:     CC-BY-NC-ND
// Created:     10/06/2015
///////////////////////////////////////

#ifndef DENSITYS2_H
#define DENSITYS2_H

#include "DataS2.h"
#include <string>
#include <vector>

class DensityS2
{
    public:
        DensityS2();
        std::vector<DataS2> get_data() const;

        double density_approx(double) const;
        double density_approx(double, unsigned int, unsigned int) const;

    protected:
        std::vector<DataS2> dataS2_load(std::string path="dataS2.txt") const;
        std::vector<DataS2> m_data;

};

#endif // DENSITYS2_H
