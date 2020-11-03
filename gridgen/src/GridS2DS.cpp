///////////////////////////////////////
// Name:        GridS2DS.cpp
// Author:      Andrzej Pisarski
// Copyright:   Andrzej Pisarski
// License:     CC-BY-NC-ND
// Created:     24/11/2019
// Modification:03/05/2020 A.P.
///////////////////////////////////////

#include "GridS2DS.h"
#include "num.h"
#include <stdexcept>
#include <cmath>
#include <algorithm>    // transform
#include <functional>   // plus, bind1st
#include <iterator>     // inserter
#include <fstream>
#include <iomanip>
#include <iostream>

#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

GridS2DS::GridS2DS(FisherRMDS const *const fmds) : m_fmds(fmds)
{
    if(m_fmds==nullptr)
    {
        std::string error="Fisher Matrix is missing.\n";
        throw std::runtime_error(error);
    }

    if(m_fmds->dim()==2 )
    {
        m_dim = 2.0;

        std::vector<std::vector<double> > a2 = a2star_pp();

        /// a2* normalized vectors (perpendicular to axis Omega_{0}^{'})
        m_a1n = num::norm(a2[0]);
        m_a2n = num::norm(a2[1]);
    }
    else{
        std::string error="Dimension of Fisher Reduced Matrix for directed searches (Grid S2) must be = 2.\n";
        throw std::domain_error(error);
    }
}

/// A2* basis perpendicular to axis Omega_{0}^{'}
std::vector<std::vector<double> > GridS2DS::a2star_pp() const
{
    ///a(n) := o(n) - hs                                // vectors perpendicular to vector 'hs'
    std::vector<double> a1 = {0., -3/2.0};
    std::vector<double> a2 = {0.,  3/2.0};

    std::vector<std::vector<double> > collect;
    collect.push_back(a1); collect.push_back(a2);

    return collect;
}

std::vector<std::vector<double> > GridS2DS::cell_vectors(double alpha, double c0=0.75,
unsigned int nfft=1048576, unsigned int data_length=344656, int mode=1) const
{
    if(data_length<=0)
    {
        std::string error="Data length can not be <= 0.\n";
        throw std::domain_error(error);
    }

    double dwp2;
    if(data_length <= nfft)
    {
        double dwp=num::delta_omega_zero_prim(c0, nfft, data_length);
        dwp2=dwp/2.0;
    }
    else{
        std::string error="Length of data need to be equal\n"
        "or less than data length of Fourier transform .\n";
        throw std::domain_error(error);
    }

    // edge's vectors:
    std::vector<double> vertex01(2);    std::vector<double> vertex02(2);

    vertex01[0]=vertex02[0]=dwp2;

    double dwp2_tan=dwp2*tan(alpha);

    for(unsigned int i=1; i<2; i++)
    {
        vertex01[i]=m_a1n[i]*dwp2_tan;
        vertex02[i]=m_a2n[i]*dwp2_tan;
    }

    std::vector<std::vector<double> > vertices;

    vertices.push_back(vertex01); vertices.push_back(vertex02);

    if(mode)
    {
        // origin:
        std::vector<double> vertex00(2);

        // vectors of remains vertex base cell
        std::vector<double> vertex03(4);

        // std::back_inserter(vertex05)
        std::transform(vertex01.begin(), vertex01.end(), vertex02.begin(), vertex03.begin(), std::plus<double>());

        vertices.push_back(vertex03);
        vertices.push_back(vertex00);
    }

    return vertices;
}


double GridS2DS::density(std::vector<double>& basis_vectors) const
{
    return num::thickness(basis_vectors, 1.0, m_dim);
    /// return std::abs( pow(M_PI, 2.0)/( 2.0*num::det(basis_vectors, 4) ) );
}

double GridS2DS::coordinate_omega1(int m, double c0, unsigned int nfft, unsigned int data_length) const
{
    double pi2 = M_PI * M_PI;
    double pi4 = pi2 * pi2;

    unsigned int k =  num::k(nfft, data_length);
    unsigned int k2 = k * k;
    unsigned int k4 = k2 * k2;

    unsigned int n = 2 * k + m;
    unsigned int n2 = n * n;
    unsigned int n4 = n2 * n2;
    unsigned int n6 = n4 * n2;

    double a = num::a(c0);
    double a2 = a * a;
    double a4 = a2 * a2;

    double e01 = - 8.0 * n6;
    double e02 = 2.0 * a2 * k2 * n2 * (-3.0 + 2.0 * n2 + n4) * pi2 ;
    double e03 = - a4 * k4 * pow(-1.0 + n2, 2) * pi4;
    double e04a = 2.0 * sqrt(n4 - a2 * k2 * (-1.0 + n2) * pi2);
    double e04b = fabs(-4.0 * n4 + a2 * k2 * (-1.0 + n4) * pi2);
    double e04 = e04a * e04b;

    double numerator = n * sqrt(e01 + e02 + e03 + e04);
    double denominator = a * k * pow(-1.0 + n2, 2) * M_PI;

    return numerator / denominator;
}

std::vector<double> GridS2DS::grid_prim(double c0, unsigned int nfft, unsigned int data_length, int mode = 1) const
{
    double DeltaOmegaZeroPrim = num::delta_omega_zero_prim(c0, nfft, data_length);
    double DeltaOmegaZeroPrim_2 = DeltaOmegaZeroPrim/2.0;

    //double DeltaOmegaZeroPrim_Min = 1.00;
    double DeltaOmegaZeroPrim_Max = 3.74;//3.74
    double Radius = 1.0;

    std::vector<double> base_prim_1vector(4, 0.0);

    base_prim_1vector[0]=DeltaOmegaZeroPrim;   //dwp2
    base_prim_1vector[1]=0.0;

    if(mode == 1 && DeltaOmegaZeroPrim <= DeltaOmegaZeroPrim_Max){
        if(DeltaOmegaZeroPrim >= 2 * Radius){
            base_prim_1vector[2] = DeltaOmegaZeroPrim_2;
            base_prim_1vector[3] = sqrt((4.0/DeltaOmegaZeroPrim)-1.0) * DeltaOmegaZeroPrim_2;
        }
        else{
            base_prim_1vector[2] = DeltaOmegaZeroPrim_2;
            base_prim_1vector[3] = sqrt(1.0 - DeltaOmegaZeroPrim_2 * DeltaOmegaZeroPrim_2) + 1.0;
        }

    }

    if(mode == 2){// Please use this algorithm only for case: 2 * data_length / nfft = k = 1, 2, 4, 8, ... .
        double thickness_cube = M_PI / 2.0;
        double thickness_hex = 2.0 * M_PI /(3.0 * sqrt(3.0));
        double a = num::a(c0);
        double a2 = a * a;
        double pi2 = M_PI * M_PI;
        double thickness_grid_old, thickness_grid = 2 * thickness_cube;

        unsigned int k = num::k(nfft, data_length);
        //std::cout << "k= " << k << std::endl;
        unsigned int k2 = k * k;

        //int m_min = 1;
        //int m_max = 9;

        std::vector<double> base_prim_test(4);

        base_prim_test[0]=a * k * M_PI;
        base_prim_test[1]=0;

        if(k == 1){
            base_prim_1vector[2] = a * M_PI / 2.0;
            base_prim_1vector[3] = 1 + sqrt(1.0 - a2 * pi2 / 4.0);
        }

        if(k >= 2 ){
            //for(int m = m_min; m < m_max; m+=2){
                unsigned int m = 1;//
                unsigned int n = 2 * k + m;
                double x = coordinate_omega1(m, c0, nfft, data_length);
                double y = (a * k * M_PI - sqrt(a2 * k2 * pi2 - 4.0 * pow(x, 2.0)))/2.0;

                base_prim_test[2] = 2.0 * (a * k * M_PI - y) / n;
                base_prim_test[3] = 2.0 * x / static_cast<double>(n);
                thickness_grid_old = thickness_grid;
                thickness_grid = density(base_prim_test);

                /*
                double thickness_grid_old = thickness_grid;
                double thickness_grid = density(base_prim_test);

                if(thickness_grid < 2 * thickness_cube && thickness_grid >= thickness_hex){
                    base_prim_1vector[2] = base_prim_test[2];
                    base_prim_1vector[3] = base_prim_test[3];
                }
                */
                if(thickness_grid < thickness_grid_old && thickness_grid >= thickness_hex){
                    base_prim_1vector[2] = base_prim_test[2];
                    base_prim_1vector[3] = base_prim_test[3];
                }
                //std::cout << "c0= " << c0 << ", omega1= " << base_prim_test[3] << ", omega0= " << base_prim_test[2]
                //<< ", DeltaOmegaZeroPrim= " << base_prim_test[0] << ", data_length= " << data_length << " , k= " << k  << ", n= " << n  << std::endl;
            //}
        }
    }

    return base_prim_1vector;
}

std::vector<double> GridS2DS::grid(double c0, double xi, unsigned int nfft, unsigned int data_length, int mode = 1) const
{
    unsigned int dim2 = m_dim * m_dim;
    std::vector<double> basis(dim2, 0.0);
    double r_scale=sqrt(1.0-c0);

    std::vector<double> base_prim_1vector=grid_prim(c0, nfft, data_length, mode);

    std::vector<double> temp;
    temp = num::multiply_AB(base_prim_1vector, num::inverse( num::cholesky( m_fmds -> postrmf(xi), m_dim), m_dim), m_dim, m_dim, m_dim);
    for(unsigned int i=0; i<dim2; i++)
        basis[i]=r_scale*temp[i];

    return basis;
}

/// Convert vector from hyper-sphere space to hyper-ellipsoid space
std::vector<double> GridS2DS::convert(double c0, double xi, const std::vector<double>& grid_prim) const
{
    unsigned int dim2 = m_dim * m_dim;
    std::vector<double> ws1(dim2, 0.0);
    double r_scale=sqrt(1.0-c0);

    std::vector<double> temp(dim2, 0.0);
    temp = num::multiply_AB(grid_prim, num::inverse( num::cholesky( m_fmds -> postrmf(xi), m_dim), m_dim), m_dim, m_dim, m_dim);
    for(unsigned int i=0; i<dim2; i++)
        ws1[i]=r_scale*temp[i];

    return ws1;
}





