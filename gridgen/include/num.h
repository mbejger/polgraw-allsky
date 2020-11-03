///////////////////////////////////////
// Name:        num.h
// Author:      Andrzej Pisarski
// Copyright:   Andrzej Pisarski
// License:     CC-BY-NC-ND
// Created:     10/06/2015
// Modification:03/05/2020
///////////////////////////////////////

#ifndef NUM_H_INCLUDED
#define NUM_H_INCLUDED

#include <vector>

namespace num
{
    static const unsigned int max_spindownDS = 3;

    double simpson(const std::vector<double>&);
    double simpson(const std::vector<double>&, const std::vector<double>&);
    double boole(const std::vector<double>&);
    double boole(const std::vector<double>&, const std::vector<double>&);

    double avg(const std::vector<double>&);
    double var(const std::vector<double>&, double avg = 0, bool unbiased = true);

    bool IsPositive (double);
    double length(const std::vector<double>&);              // return length of the a vector
    std::vector<double> norm(const std::vector<double>&);   // return normalized a vector

    std::vector<double> minor1(const std::vector<double>&, unsigned int, unsigned int, unsigned int);
    double det(const std::vector<double>, unsigned int);

    std::vector<double> cofactor(const std::vector<double>&, unsigned int);
    std::vector<double> transpose(const std::vector<double>&, unsigned int, unsigned int);
    std::vector<double> inverse(const std::vector<double>&, unsigned int);

    std::vector<double> cholesky(const std::vector<double>&, unsigned int);
    std::vector<double> multiply_A(const std::vector<double>&, double);//multiply a matrix A (vector) by a constant
    std::vector<double> multiply_AB(const std::vector<double>&, const std::vector<double>&,
                             unsigned int, unsigned int, unsigned int);//multiply a matrix A by a matrix B
    std::vector<double> diagonal(double, unsigned int);

    int sign(unsigned int, unsigned int);                   // 1 if sum is even, -1 if is not
    int sign(double);                                       // return sign of double
    double epsilon();                                       // precision == 10^-12
    double chop(double, double precision=num::epsilon());
    void chop(std::vector<double>&, double precision=num::epsilon());

    /// Matrix rotation functions:
    std::vector<double> scb(const std::vector<double>&);
    std::vector<double> rotation(const std::vector<double>&, int, int, unsigned int);  //rotation right
    std::vector<double> single_rot(const std::vector<double>&, int, int, int, unsigned int);
    std::vector<double> rot(const std::vector<double>&, unsigned int dim=4);//makes multi rotations
    // (using 'single_rot' T_{dim-1} times [T - triangular number]) to get lower triangular matrix

    /// Fourier frequency (_prim - space with hyper-spheres)
    double delta_omega_zero(unsigned int, unsigned int);
    double delta_omega_zero_prim(double, unsigned int, unsigned int);
    double a(double c0);    //Coefficient 'a' - can be used to translate from original space (ellipses)
    unsigned int k(unsigned int, unsigned int); // for more info see 'num.cpp'

    /// Sphere packings:
    std::vector<double> gram(unsigned int);//return a Gram matrix (for a given dimension);
    std::vector<double> a_star(unsigned int);//return A[n]* matrix (n == dimension);
    std::vector<double> a_star_r1(unsigned int);//return A[n]* matrix; radius of covering is equal to 1;

    double sphere_volume(unsigned int dim, double radius=1.0);
    double covering_radius_a_star(unsigned int);//covering radius for a grid a_star;
    double thickness(const std::vector<double>&, double, unsigned int);

}

#endif // NUM_H_INCLUDED
