///////////////////////////////////////
// Name:        FisherRM.h
// Author:      Andrzej Pisarski
// Copyright:   Andrzej Pisarski
// License:     CC-BY-NC-ND
// Created:     20/05/2015
// Modification:28/08/2017 A.P.
///////////////////////////////////////

#ifndef FISHERRM_H
#define FISHERRM_H

#include <string>
#include <vector>

class FisherRM
{
    public:
        FisherRM();
        FisherRM(const std::vector<double>& ephemeris1, const std::vector<double>& ephemeris2, int spindown_=1);
        // vector<double>& ephemeris1, ephemeris2

        std::vector<double> get_prefisher() const;
        std::vector<double> postrmf(double) const;
        unsigned int get_ephemeris_length() const;
        unsigned int dim() const;                   // gives dimension of Fisher matrix

    protected:
        std::vector<double> m_prefisher;            // (unchanged) elements of reduced Fisher matrix
        const unsigned int m_ephemeris_length;      // number of elements (data length) in ephemeris1
                                                    // (ephemeris2 should have this same length)
        int spindown;
};

#endif // FISHERRM_H
