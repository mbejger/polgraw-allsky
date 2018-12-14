///////////////////////////////////////
// Name:        ReadEphemeris.h
// Author:      Andrzej Pisarski
// Copyright:   Andrzej Pisarski
// License:     CC-BY-NC-ND
// Created:     20/05/2015
///////////////////////////////////////

#ifndef READEPHEMERIS_H
#define READEPHEMERIS_H

#include <vector>
#include <string>


class ReadEphemeris
{
    public:
        ReadEphemeris(const std::vector<std::string>);     /// DEFAULT C-TOR ?? !!!

        const std::vector<double>& get_ephemeris1() const; // \mu_{1} - return Cartesian coordinates at SSB
        const std::vector<double>& get_ephemeris2() const; // \mu_{2} - return Cartesian coordinates at SSB
        unsigned int get_length() const;

    private:
        std::vector<double> m_ephemeris1;         // \mu_{1} - Cartesian coordinates at SSB.
        std::vector<double> m_ephemeris2;         // \mu_{2} - Cartesian coordinates at SSB.
        unsigned int m_length;                    // number of elements in ephemeris1, ephemeris2

        std::vector<char> m_flag; // files streams states
        double m_epsm;            // angle of ecliptic inclination
        double m_ce;              // cosine of an angle of ecliptic inclination
        std::vector<double> m_DetSSB;
};

#endif // READEPHEMERIS_H
