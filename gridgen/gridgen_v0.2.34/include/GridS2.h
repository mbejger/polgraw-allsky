///////////////////////////////////////
// Name:        GridS2.h
// Author:      Andrzej Pisarski
// Copyright:   Andrzej Pisarski
// License:     CC-BY-NC-ND
// Created:     10/06/2015
///////////////////////////////////////

#ifndef GRIDS2_H
#define GRIDS2_H

#include "FisherRM.h"
#include "CellS2.h"
#include "DataS2.h"
#include <string>
#include <vector>

class GridS2
{
    public:
        GridS2(FisherRM const *const, int, int);

        std::vector<double> grid_prim(double, unsigned int, bool, bool, std::string path="dataS2.txt") const;
        std::vector<double> grid(double, double, unsigned int, bool, bool, std::string path="dataS2.txt") const;
        std::vector<double> convert(double, double, const std::vector<double>&) const; // convert vector
                                                // from hyper-sphere space to hyper-ellipsoid space
        std::vector<std::vector<double> > a4star_pp() const;    // A4* basis (vector of 4 base vectors
                                                                // perpendicular to axis Omega_{0}^{'})
        std::vector<std::vector<double> > cell_vectors(double, double, unsigned int, int) const;
                                            // 4-dimensional
                                            // (int = 0): 4 cell vectors (4 edges);
                                            // (int = 1 - default set): 16 cell vectors
                                            // (+ 1 origin = {0,0,0,0} - last vector))
        //double density(CellS2& );
        double density(std::vector<double>& ) const; // need 4 base vectors in 4D space with hyper-spheres
        CellS2 find_alpha(double, double, double, unsigned int) const;

    protected:
        FisherRM const *const m_fm;
        const int m_NAlpha;                 // incrementation deep of root finding algorithm
        const int m_NRadius;                // incrementation deep of covering radius (deep hole)
                                            // founding algorithm


        /// a4* normalized vectors
        std::vector<double> m_a1n;
        std::vector<double> m_a2n;
        std::vector<double> m_a3n;
        std::vector<double> m_a4n;

        std::vector<DataS2> dataS2_load(std::string path="dataS2.txt") const;
        void                dataS2_save(std::vector<DataS2>&, std::string path) const;

    private:
        CellS2 new_radius_approx(double, std::vector<double>&, std::vector<std::vector<double> >&,
                                 std::vector<double> &) const;
        CellS2 find_radius(double, double, unsigned int) const;
        int check_nalpha(int);
        int check_nradius(int);

};

#endif // GRIDS2_H
