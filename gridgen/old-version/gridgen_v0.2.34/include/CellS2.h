///////////////////////////////////////
// Name:        CellS2.h
// Author:      Andrzej Pisarski
// Copyright:   Andrzej Pisarski
// License:     CC-BY-NC-ND
// Created:     10/06/2015
///////////////////////////////////////

#ifndef CELLS2_H
#define CELLS2_H

#include <vector>

class CellS2
{
    public:
        CellS2();
        CellS2(double, double, std::vector<double>);
        //CellS2(const CellS2&);

        //void operator=(const CellS2&);

        double m_alpha;                   // angle between omega0prim axis and (any)one of basis vectors
        double m_radius;                  // radius of covering
        std::vector<double> get_hole() const;   // return m_hole

    private:
        std::vector<double> m_hole;     // vector between origin and (deep) hole approximation
};

bool compare_CellS2_alpha(const CellS2&, const CellS2&);
bool compare_CellS2_radius(const CellS2&, const CellS2&);
bool compare_CellS2_alpha_decrease(const CellS2&, const CellS2&);
bool compare_CellS2_radius_decrease(const CellS2&, const CellS2&);

#endif // CELLS2_H
