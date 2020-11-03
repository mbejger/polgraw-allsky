///////////////////////////////////////
// Name:        FisherRMDS.h
// Author:      Andrzej Pisarski
// Copyright:   Andrzej Pisarski
// License:     CC-BY-NC-ND
// Created:     13/10/2019
// Modification:08/12/2019 A.P.
///////////////////////////////////////

#ifndef FISHERRMDS_H
#define FISHERRMDS_H

#include <string>
#include <vector>

class FisherRMDS
{
    public:
        FisherRMDS();
        FisherRMDS(int spindown_=1);

        std::vector<double> postrmf(double) const;
        unsigned int dim() const;                   // gives dimension of Fisher matrix (directed searches case)
        unsigned int get_spindown() const;

    protected:
        std::vector<double> fisherMDS(double, int) const; // like 'postrmf' but takes also spindown value (int)
        std::vector<double> fisher_elemetsDS(double, unsigned int) const; /// Auxiliary function - simple way to add Fisher matrix elements (in directed searches case).
        unsigned int m_spindown;
};

#endif // FISHERRMDS_H
