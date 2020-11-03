///////////////////////////////////////
// Name:        GridS1DS.cpp
// Author:      Andrzej Pisarski
// Copyright:   Andrzej Pisarski
// License:     CC-BY-NC-ND
// Created:     25/10/2019
// Modification:08/12/2019 A.P.
///////////////////////////////////////

#include "GridS1DS.h"
#include "num.h"
#include <stdexcept>
#include <cmath>
#include <string>

#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

GridS1DS::GridS1DS(const FisherRMDS *const fmds, unsigned int spindown_): DensityS1DS(spindown_), m_fmds(fmds)
{
    if(m_fmds==nullptr)
    {
        std::string error="Fisher Matrix for Directed Searches is missing.\n";
        throw std::runtime_error(error);
    }
}

/// Produce grid in hyper-ellipsoid space;
/// be careful and do not convert vectors twice (exp. by using 'convert' function)
std::vector<double> GridS1DS::grid(double c0, double xi, unsigned int nfft, unsigned int data_length) const
{
    unsigned int dim2 = m_dim * m_dim;
    std::vector<double> ws1(dim2, 0.0);
    double r_scale=sqrt(1.0-c0);

    std::vector<double> q0 = DensityS1DS::grid_prim(c0, nfft, data_length);
    //ws1[0]=DensityS1::density(q0);
    std::vector<double> temp;
    temp = num::multiply_AB(q0, num::inverse( num::cholesky( m_fmds->postrmf(xi), m_dim), m_dim), m_dim, m_dim, m_dim);
    for(unsigned int i=0; i<dim2; i++)
        ws1[i]=r_scale*temp[i]; //s1[i+1]=

    ///ws1 = DensityS1DS::grid_prim(c0, nfft, data_length);

    return ws1;
}

// convert vector from hyper-sphere space to hyper-ellipsoid space
std::vector<double> GridS1DS::convert(double c0, double xi, const std::vector<double>& grid_prim)
{
    unsigned int dim2 = m_dim * m_dim;
    std::vector<double> ws1(dim2);
    double r_scale=sqrt(1.0-c0);

    std::vector<double> temp(dim2, 0.0);
    temp = num::multiply_AB(grid_prim, num::inverse( num::cholesky( m_fmds->postrmf(xi), m_dim), m_dim), m_dim, m_dim, m_dim);

    for(unsigned int i=0; i<dim2; i++)
        ws1[i]=r_scale*temp[i];
    //test:
    //return num::cholesky( m_fmds->postrmf(xi), m_dim);
    return ws1;
}

/// double GridS1DS::density(double c0,  unsigned int nfft, unsigned int data_length) const
/// {
///    return DensityS1DS::density(c0, nfft, data_length);
/// }


