///////////////////////////////////////
// Name:        GridS2DS.h
// Author:      Andrzej Pisarski
// Copyright:   Andrzej Pisarski
// License:     CC-BY-NC-ND
// Created:     02/11/2019
// Modification:11/04/2020 A.P.
///////////////////////////////////////

#ifndef GRIDS2DS_H
#define GRIDS2DS_H

#include "FisherRMDS.h"
#include "CellS2.h"
#include "DataS2.h"
#include <string>
#include <vector>

class GridS2DS
{
    public:
        GridS2DS(FisherRMDS const * const );

        std::vector<double> grid_prim(double, unsigned int, unsigned int, int) const;
        std::vector<double> grid(double, double, unsigned int, unsigned int, int) const;
        std::vector<double> convert(double, double, const std::vector<double>&) const; // convert vector
                                                // from hyper-sphere space to hyper-ellipsoid space
        std::vector<std::vector<double> > a2star_pp() const;    // A2* basis (vector of 2 base vectors
                                                                // perpendicular to axis Omega_{0}^{'})
        std::vector<std::vector<double> > cell_vectors(double, double, unsigned int, unsigned int, int) const;
                                            // 2-dimensional
                                            // (int = 0): 2 cell vectors (2 edges);
                                            // (int = 1 - default set): 16 cell vectors
                                            // (+ 1 origin = {0,0,0,0} - last vector))
        //double density(CellS2& );
        double density(std::vector<double>& ) const; // need 2 base vectors in 2D space with hyper-spheres

    protected:
        FisherRMDS const * const m_fmds;
        unsigned int m_dim;                 // An* dimension

        /// a2* normalized vectors
        std::vector<double> m_a1n;
        std::vector<double> m_a2n;

        double coordinate_omega1(int, double, unsigned int, unsigned int) const;
};

#endif // GRIDS2DS_H
