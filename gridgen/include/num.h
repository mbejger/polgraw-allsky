///////////////////////////////////////
// Name:        num.h
// Author:      Andrzej Pisarski
// Copyright:   Andrzej Pisarski
// License:     CC-BY-NC-ND
// Created:     10/06/2015
///////////////////////////////////////

#ifndef NUM_H_INCLUDED
#define NUM_H_INCLUDED

#include <vector>

namespace num
{
    double simpson(const std::vector<double>&);
    double simpson(const std::vector<double>&, const std::vector<double>&);

    bool IsBTZ (double);
    double length(const std::vector<double>&);              // return length of the a vector
    std::vector<double> norm(const std::vector<double>&);   // return normalized a vector

    std::vector<double> minor2(const std::vector<double>&, unsigned int, unsigned int, unsigned int);
    double det(const std::vector<double>, unsigned int);

    std::vector<double> cofactor(const std::vector<double>&, unsigned int);
    std::vector<double> transpose(const std::vector<double>&, unsigned int, unsigned int);
    std::vector<double> inverse(const std::vector<double>&, unsigned int);

    std::vector<double> cholesky(const std::vector<double>&, unsigned int);
    std::vector<double> mult(const std::vector<double>&, const std::vector<double>&,
                             unsigned int, unsigned int, unsigned int);
    std::vector<double> diagonal(double, unsigned int);

    int sign(unsigned int, unsigned int);                   // 1 if sum is even, -1 if is not
    int sign(double);                                       // return sign of double
    double epsilon();                                       // precision == 10^-12
    double chop(double, double precision=num::epsilon());
    void chop(std::vector<double>&, double precision=num::epsilon());

    /// Matrix rotation functions:
    std::vector<double> scb(const std::vector<double>&);
    std::vector<double> rotright(const std::vector<double>&, int);   //rotation right
    std::vector<double> rotleft(const std::vector<double>&, int);    //rotation left
    std::vector<double> rot(const std::vector<double>&, int, int, int);

    /// Fourier frequency (_prim - space with hyper-spheres)
    double delta_omega_zero(unsigned int, unsigned int);
    double delta_omega_zero_prim(double, unsigned int, unsigned int);

    /// Experimental functions:
    double simpson_adv(const std::vector<double>&);
    double simpson_adv(const std::vector<double>&, const std::vector<double>&);
    double sum(std::vector<double>&);
}

#endif // NUM_H_INCLUDED

/// There is no need to set function 'const' which have global scope
/// Need to use only with member functions
/// A const variable has to be declared within the class, but it cannot be defined in it.
