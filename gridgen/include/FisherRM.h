///////////////////////////////////////
// Name:        FisherRM.h
// Author:      Andrzej Pisarski
// Copyright:   Andrzej Pisarski
// License:     CC-BY-NC-ND
// Created:     20/05/2015
// Modification:31/05/2020 A.P.
///////////////////////////////////////

#ifndef FISHERRM_H
#define FISHERRM_H

#include <string>
#include <vector>

class FisherRM
{
    public:
        FisherRM();
        // vector<double>& ephemeris1, ephemeris2
        FisherRM(const std::vector<double>& ephemeris1, const std::vector<double>& ephemeris2, int spindown_=1);
        FisherRM(const std::vector<std::vector<double> >& ephemeris, const std::vector<double>& sigma, int spindown_=1);

        std::vector<double> get_prefisher() const;
        std::vector<double> get_prefisher(int) const;
        std::vector<double> postrmf(double) const;
        std::vector<double> postrmf(double, int) const;
        std::vector<double> get_sigma() const;
        unsigned int get_ephemeris_length() const;
        unsigned int dim() const;                   // gives dimension of Fisher matrix

    protected:
        //std::vector<double> m_prefisher;     // @depreciated
        std::vector<std::vector<double> > m_prefisher_net;// vector of reduced Fisher matrix(s) (invariable elements)
        const unsigned int m_ephemeris_length;      // number of elements (data length) in ephemeris[0]
                                                    // (ephemeris[0] == ephemeris1, ephemeris[1] == ephemeris2, ... )
                                                    // (ephemeris[x] should have this same length as ephemeris1).
        std::vector<double> m_sigma;
        const unsigned int m_number_of_detectors;
        int m_spindown;

};

#endif // FISHERRM_H
