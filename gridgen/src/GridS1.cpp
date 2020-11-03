///////////////////////////////////////
// Name:        GridS1.cpp
// Author:      Andrzej Pisarski
// Copyright:   Andrzej Pisarski
// License:     CC-BY-NC-ND
// Created:     10/06/2015
// Modification:07/09/2017 A.P.
///////////////////////////////////////

#include "GridS1.h"
#include "num.h"
#include <stdexcept>
#include <cmath>
#include <string>

#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

GridS1::GridS1(const FisherRM *const fm, unsigned int spindown_): DensityS1(spindown_), m_fm(fm)
{
    if(m_fm==nullptr)
    {
        std::string error="Fisher Matrix is missing.\n";
        throw std::runtime_error(error);
    }
}

/// Produce grid in hyper-ellipsoid space;
/// be careful and do not convert vectors twice (exp. by using 'convert' function)
std::vector<double> GridS1::grid(double c0, double xi, unsigned int nfft) const
{
    unsigned int dim2 = m_dim * m_dim;
    std::vector<double> ws1(dim2, 0.0);
    unsigned int data_length = m_fm->get_ephemeris_length();

    double r_scale=sqrt(1.0-c0);
    std::vector<double> q0 = DensityS1::grid_prim(c0, nfft, data_length);
    //ws1[0]=DensityS1::density(q0);
    std::vector<double> temp;
    temp = num::multiply_AB(q0, num::inverse( num::cholesky( m_fm->postrmf(xi), m_dim), m_dim), m_dim, m_dim, m_dim);
    for(unsigned int i=0; i<dim2; i++)
        ws1[i]=r_scale*temp[i]; //s1[i+1]=

    ///ws1 = DensityS1::grid_prim(c0, nfft, data_length);

    return ws1;
}

// convert vector from hyper-sphere space to hyper-ellipsoid space
std::vector<double> GridS1::convert(double c0, double xi, const std::vector<double>& grid_prim)
{
    int dim2 = m_dim * m_dim;
    std::vector<double> ws1(dim2);
    double r_scale=sqrt(1.0-c0);

    std::vector<double> temp(dim2, 0.0);
    temp = num::multiply_AB(grid_prim, num::inverse( num::cholesky( m_fm->postrmf(xi), m_dim), m_dim), m_dim, m_dim, m_dim);

    for(int i=0; i<dim2; i++)
        ws1[i]=r_scale*temp[i]; //ws1[i]=

    return ws1;
}

double GridS1::density(double c0,  unsigned int nfft) const
{
    unsigned int data_length=m_fm->get_ephemeris_length();
    return DensityS1::density(c0, nfft, data_length);
}

///   double GridS1::density(double c0,  unsigned int nfft, unsigned int data_length) const
///   {
///        return DensityS1::density(c0, nfft, data_length);
///   }

