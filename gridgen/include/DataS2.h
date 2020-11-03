///////////////////////////////////////
// Name:        DataS2.h
// Author:      Andrzej Pisarski
// Copyright:   Andrzej Pisarski
// License:     CC-BY-NC-ND
// Created:     16/06/2015
///////////////////////////////////////

#ifndef DATAS2_H
#define DATAS2_H

/// Structure 'DataS2' keep data: (delta_omega0_prim, density, alpha, n_root, n_sampling)
class DataS2
{
    public:
        DataS2();
        DataS2(double, double, double, double, int, int);

        double m_omega0p;
        double m_density;
        double m_alpha;
        double m_radius;

        int get_nalpha() const;     // return incrementation deep (m_NAlpha) of root finding algorithm
        int get_nradius() const;    // return incrementation deep (m_NRadius)
                                    // of deep hole founding algorithm

    private:
        int m_NAlpha;
        int m_NRadius;

};

bool compare_DataS2_omega0p(const DataS2&, const DataS2&);
bool compare_DataS2_density(const DataS2&, const DataS2&);
bool compare_DataS2_alpha(const DataS2&, const DataS2&);
bool compare_DataS2_radius(const DataS2&, const DataS2&);

///decrease
bool compare_DataS2_omega0p_decrease(const DataS2&, const DataS2&);
bool compare_DataS2_density_decrease(const DataS2&, const DataS2&);
bool compare_DataS2_alpha_decrease(const DataS2&, const DataS2&);
bool compare_DataS2_radius_decrease(const DataS2&, const DataS2&);

#endif // DATAS2_H
