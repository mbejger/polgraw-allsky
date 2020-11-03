///////////////////////////////////////
// Name:        GridS2.cpp
// Author:      Andrzej Pisarski
// Copyright:   Andrzej Pisarski
// License:     CC-BY-NC-ND
// Created:     10/06/2015
///////////////////////////////////////

#include "GridS2.h"
#include "DataS2.h"
#include "num.h"
#include <stdexcept>
#include <cmath>
#include <algorithm>    // transform
#include <functional>   // plus, bind1st
#include <iterator>     // inserter
#include <fstream>
#include <iomanip>
//#include <iostream>

#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

GridS2::GridS2(FisherRM const *const fm, int n_alpha=35, int n_radius=20):
    m_fm(fm), m_NAlpha( check_nalpha(n_alpha) ), m_NRadius( check_nradius(n_radius) )
{
    if(m_fm==nullptr)
    {
        std::string error="Fisher Matrix is missing.\n";
        throw std::runtime_error(error);
    }

    if(m_fm->dim()==4)/// || m_fm->dim()==5
    {
        m_dim = 4;
        std::vector<std::vector<double> > a4 = a4star_pp();

        /// a4* normalized vectors (perpendicular to axis Omega_{0}^{'})
        m_a1n = num::norm(a4[0]);
        m_a2n = num::norm(a4[1]);
        m_a3n = num::norm(a4[2]);
        m_a4n = num::norm(a4[3]);
    }
    else{
        std::string error="Dimension of Fisher Reduced Matrix must be = 4.\n";
        throw std::domain_error(error);
    }
}

/// A4* basis perpendicular to axis Omega_{0}^{'}
std::vector<std::vector<double> > GridS2::a4star_pp() const
{
    //double sq2=sqrt(2);
    //double sq2_d=2.0*sq2;
    //double sq3=sqrt(3);
    double sq5=sqrt(5);
    double sq5_3=sqrt(5.0/3.0);
    double sq5_6=sqrt(5.0/6.0);

    //std::vector<double> o1 = {1.0/sq2_d, -sq5_6/2.0, -sq5_3/2.0, -sq5/2.0};
    //std::vector<double> o2 = {1.0/sq2_d, sqrt(7.5)/2.0, 0., 0. };   //sqrt(15.0/2.0)
    //std::vector<double> o3 = {1.0/sq2_d, -sq5_6/2.0, sq5_3, 0.};
    //std::vector<double> o4 = {1.0/sq2_d, -sq5_6/2.0, -sq5_3/2.0, sq5/2.0};
    //std::vector<double> hs = {1.0/sq2_d, 0., 0., 0.};   // symmetry axis for vector o(n) (o1+o2+o3+o4)/4
    //std::vector<double> hsn = {1.0, 0., 0., 0.};        // normalized 'hs' vector

    ///a(n) := o(n) - hs                                // vectors perpendicular to vector 'hs'
    std::vector<double> a1 = {0., -sq5_6/2.0, -sq5_3/2.0, -sq5/2.0};
    std::vector<double> a2 = {0., sqrt(7.5)/2.0, 0., 0. };     // sqrt(15.0/2.0)
    std::vector<double> a3 = {0., -sq5_6/2.0, sq5_3, 0.};
    std::vector<double> a4 = {0., -sq5_6/2.0, -sq5_3/2.0, sq5/2.0};

    std::vector<std::vector<double> > collect;
    collect.push_back(a1); collect.push_back(a2);
    collect.push_back(a3); collect.push_back(a4);

    return collect;
}

std::vector<std::vector<double> > GridS2::cell_vectors(double alpha, double c0=0.75,
                                                       unsigned int nfft=1048576, int mode=1) const
{
    unsigned int data_length = m_fm->get_ephemeris_length();
    if(data_length<=0)
    {
        std::string error="Data length can not be <= 0.\n";
        throw std::domain_error(error);
    }

    double dwp4;
    if(data_length <= nfft)
    {
        double dwp=num::delta_omega_zero_prim(c0, nfft, data_length);
        dwp4=dwp/4.0;
    }
    else{
        std::string error="Length of Fourier transform need to be equal\n"
        "or longer than data length.\n";
        throw std::domain_error(error);
    }

    // edge's vectors:
    std::vector<double> vertex01(4);    std::vector<double> vertex02(4);
    std::vector<double> vertex03(4);    std::vector<double> vertex04(4);

    vertex01[0]=vertex02[0]=vertex03[0]=vertex04[0]=dwp4;

    double dwp4_tan=dwp4*tan(alpha);

    for(unsigned int i=1; i<4; i++)
    {
        vertex01[i]=m_a1n[i]*dwp4_tan;
        vertex02[i]=m_a2n[i]*dwp4_tan;
        vertex03[i]=m_a3n[i]*dwp4_tan;
        vertex04[i]=m_a4n[i]*dwp4_tan;
    }

    std::vector<std::vector<double> > vertices;

    vertices.push_back(vertex01); vertices.push_back(vertex02);
    vertices.push_back(vertex03); vertices.push_back(vertex04);

    if(mode)
    {
        // origin:
        std::vector<double> vertex00(4);

        // vectors of remains vertex base cell
        std::vector<double> vertex05(4);  std::vector<double> vertex06(4);
        std::vector<double> vertex07(4);  std::vector<double> vertex08(4);
        std::vector<double> vertex09(4);  std::vector<double> vertex10(4);
        std::vector<double> vertex11(4);  std::vector<double> vertex12(4);
        std::vector<double> vertex13(4);  std::vector<double> vertex14(4);
        std::vector<double> vertex15(4);

        // std::back_inserter(vertex05)
        std::transform(vertex01.begin(), vertex01.end(), vertex02.begin(), vertex05.begin(), std::plus<double>());
        std::transform(vertex01.begin(), vertex01.end(), vertex03.begin(), vertex06.begin(), std::plus<double>());
        std::transform(vertex01.begin(), vertex01.end(), vertex04.begin(), vertex07.begin(), std::plus<double>());
        std::transform(vertex02.begin(), vertex02.end(), vertex03.begin(), vertex08.begin(), std::plus<double>());
        std::transform(vertex02.begin(), vertex02.end(), vertex04.begin(), vertex09.begin(), std::plus<double>());
        std::transform(vertex03.begin(), vertex03.end(), vertex04.begin(), vertex10.begin(), std::plus<double>());
        std::transform(vertex01.begin(), vertex01.end(), vertex08.begin(), vertex11.begin(), std::plus<double>());
        std::transform(vertex04.begin(), vertex04.end(), vertex08.begin(), vertex12.begin(), std::plus<double>());
        std::transform(vertex03.begin(), vertex03.end(), vertex07.begin(), vertex13.begin(), std::plus<double>());
        std::transform(vertex02.begin(), vertex02.end(), vertex07.begin(), vertex14.begin(), std::plus<double>());
        std::transform(vertex07.begin(), vertex07.end(), vertex08.begin(), vertex15.begin(), std::plus<double>());

        vertices.push_back(vertex05); vertices.push_back(vertex06);
        vertices.push_back(vertex07); vertices.push_back(vertex08);
        vertices.push_back(vertex09); vertices.push_back(vertex10);
        vertices.push_back(vertex11); vertices.push_back(vertex12);
        vertices.push_back(vertex13); vertices.push_back(vertex14);
        vertices.push_back(vertex15);
        vertices.push_back(vertex00);
    }

    return vertices;
}

CellS2 GridS2::new_radius_approx(double alpha, std::vector<double>& shrink,
                                 std::vector<std::vector<double> >& vertices,
                                 std::vector<double>& hole) const
{
    std::vector<double> v1 = vertices[0]; std::vector<double> v2 = vertices[1];
    std::vector<double> v3 = vertices[2]; std::vector<double> v4 = vertices[3];

    std::transform(v1.begin(), v1.end(), v1.begin(), std::bind1st(std::multiplies<double>(),shrink[0]));
    std::transform(v2.begin(), v2.end(), v2.begin(), std::bind1st(std::multiplies<double>(),shrink[1]));
    std::transform(v3.begin(), v3.end(), v3.begin(), std::bind1st(std::multiplies<double>(),shrink[2]));
    std::transform(v4.begin(), v4.end(), v4.begin(), std::bind1st(std::multiplies<double>(),shrink[3]));

    std::transform(v1.begin(), v1.end(), v2.begin(), v1.begin(), std::plus<double>());
    std::transform(v3.begin(), v3.end(), v4.begin(), v3.begin(), std::plus<double>());
    std::transform(v1.begin(), v1.end(), v3.begin(), v1.begin(), std::plus<double>());
    std::transform(v1.begin(), v1.end(), hole.begin(), v1.begin(), std::plus<double>());

    std::vector<CellS2> search_hole;
    for(int i=0; i<16; i++)
    {

        std::vector<double> radius;
        std::transform(v1.begin(), v1.end(), vertices[i].begin(), inserter(radius, radius.begin()), std::minus<double>());
        search_hole.push_back( CellS2(alpha, num::length(radius), v1 ) );
    }

    //sort(search_hole.begin(), search_hole.end(), compare_CellS2_radius);
    //return search_hole[0];

    std::vector<CellS2>::const_iterator it;
    it = std::min_element(search_hole.begin(), search_hole.end(), compare_CellS2_radius );

    return *it;
}

CellS2 GridS2::find_radius(double alpha, double c0, unsigned int nfft) const
{
    double sa=-1./3., sb=1./3., dd=1./6.;
    std::vector<std::vector<double> > all_cell_vectors = cell_vectors(alpha, c0, nfft);
    std::vector<double> hole(4,0.);   // approximation of deep hole (start point)

    for(int i=0; i<4; i++)
        std::transform(hole.begin(), hole.end(), all_cell_vectors[i].begin(), hole.begin(), std::plus<double>());

    std::transform(hole.begin(), hole.end(), hole.begin(), std::bind1st(std::multiplies<double>(),0.5));
    //for(int i=0; i<4; i++)
    //    hole[i]*=0.5;

    CellS2 to_find;
    for(int n=0; n<=m_NRadius; n++)
    {
        for(double i=sa; i<=sb; i+=dd)
        {
            for(double j=sa; j<=sb; j+=dd)
            {
                for(double k=sa; k<=sb; k+=dd)
                {
                    for(double o=sa; o<=sb; o+=dd)
                    {
                        std::vector<double> koji={i,j,k,o};     // koji - light, happiness
                        CellS2 temp=new_radius_approx(alpha, koji, all_cell_vectors, hole);
                        if( temp.m_radius>to_find.m_radius )
                        {
                            to_find=temp;
                        }
                    }
                }
            }

        }
        hole = to_find.get_hole();
        sa*=1./3.;sb*=1./3.;dd*=1./3.;
    }

    return to_find;

}

CellS2 GridS2::find_alpha(double alpha_min, double alpha_max, double c0, unsigned int nfft) const
{
    double alpha_left=alpha_min, alpha_right=alpha_max, alpha_diff;
    CellS2 middle;
    CellS2 left;

    for(int i=0; i<m_NAlpha; i++)
    {
        alpha_diff=std::abs(alpha_right - alpha_left);
        double alpha_middle = alpha_left + 0.5*alpha_diff;

        left = find_radius(alpha_left, c0, nfft);
        double radius_left = left.m_radius;

        middle = find_radius(alpha_middle, c0, nfft);
        double radius_middle = middle.m_radius;    //+ num::epsilon: radius of covering = 1 - epsilon

        //std::cout << std::setprecision(12) << "alpha_left=" << alpha_left
            //<< " (" << alpha_left/(M_PI/2.0) << ")"
        //<< ", radius_left=" << left.m_radius <<"\n";
        //std::cout << std::setprecision(12) << "alpha_middle=" << alpha_middle
            //<< " (" << alpha_middle/(M_PI/2.0) << ")"
        //<< ", radius_midle=" << middle.m_radius <<"\n\n";

        if( num::sign(radius_left-1.)!=num::sign(radius_middle-1.) )
            alpha_right = alpha_middle;
        else
            alpha_left = alpha_middle;

    }

    if(middle.m_radius-1.<=num::epsilon())
        return middle;
    else
        return left;
}

/*
double GridS2::density(CellS2& root = find_alpha(0.001, M_PI/2.) )
{
    std::vector <double> basis_vectors = root.get_hole();
    return std::abs( pow(M_PI, 2.0)/( 2.0*num::det(basis_vectors, 4) ) );
}
*/
double GridS2::density(std::vector<double>& basis_vectors) const
{
    return num::thickness(basis_vectors, 1.0, m_dim);
    /// return std::abs( pow(M_PI, 2.0)/( 2.0*num::det(basis_vectors, 4) ) );
}

std::vector<double> GridS2::grid_prim(double c0, unsigned int nfft,
        bool s2use, bool s2save, std::string path) const
{
    double DeltaOmZeroPrim = num::delta_omega_zero_prim(c0, nfft, m_fm->get_ephemeris_length() );
    double DeltaOmZeroPrim_min = DeltaOmZeroPrim - num::epsilon();
    double DeltaOmZeroPrim_max = DeltaOmZeroPrim + num::epsilon();

    std::vector<DataS2> data_vec;
    DataS2 data;
    CellS2 root;
    bool is_found = false;
    bool ii_is_greater=false;   // ii - increment index
    unsigned int found_index=0;

    //std::cout << "c0=" << c0 << " DelOmZP=" << DeltaOmZeroPrim << "\n";
    // std::cout << "s2use=" << s2use << ", s2save=" << s2save << "\n";

    if(s2use)
    {
        data_vec = dataS2_load();
        for(unsigned int i=0; i<=data_vec.size(); i++)
            if( data_vec[i].m_omega0p<=DeltaOmZeroPrim_max && data_vec[i].m_omega0p >=DeltaOmZeroPrim_min )
            {
                data = data_vec[i];
                is_found = true;
                found_index=i;
            }

        if(is_found)
        {
            if(data.get_nalpha()<m_NAlpha || data.get_nradius()<m_NRadius)
            {
                ii_is_greater = true;
                root = find_alpha(0.001, M_PI/2., c0, nfft);
            }
        }
        else{
            root = find_alpha(0.001, M_PI/2., c0, nfft);
        }
    }
    else{
        root = find_alpha(0.001, M_PI/2., c0, nfft);
    }

    std::vector<std::vector<double> > basis_prim;
    if(is_found && !ii_is_greater)
        basis_prim = cell_vectors(data.m_alpha, c0, nfft, 0);
    else
        basis_prim = cell_vectors(root.m_alpha, c0, nfft, 0);

    std::vector<double> base_prim_1vector(16);
    for(int i=1; i<4; i++)
        for(int j=0; j<4; j++)
            base_prim_1vector[4*i+j]=basis_prim[i][j];

    //for(int i=0; i<4; i++)
    //    base_prim_1vector[0]+=basis_prim[i][0];

    base_prim_1vector[0]=DeltaOmZeroPrim;   //dwp4

    double ro = density(base_prim_1vector);
    if(s2save)
    {
        if(s2use)
        {
            if(is_found)
            {
                if(ii_is_greater) // for replace
                {
                    double alpha = root.m_alpha;
                    double radius = root.m_radius;
                    data_vec.at(found_index)=DataS2(DeltaOmZeroPrim, ro, alpha, radius, m_NAlpha, m_NRadius);
                    dataS2_save(data_vec, path);
                }
            }
            else{
                double alpha = root.m_alpha;
                double radius = root.m_radius;
                data_vec.push_back(DataS2(DeltaOmZeroPrim, ro, alpha, radius, m_NAlpha, m_NRadius));
                sort(data_vec.begin(), data_vec.end(), compare_DataS2_omega0p);
                dataS2_save(data_vec, path);
            }
        }
        else{
            std::vector<DataS2> data_vec2;
            data_vec2 = dataS2_load();
            DataS2 data2;
            bool is_found2 = false;
            unsigned int found_index2 = 0;
            for(unsigned int i=0; i<=data_vec.size(); i++)
            if( data_vec2[i].m_omega0p<=DeltaOmZeroPrim_max && data_vec2[i].m_omega0p >=DeltaOmZeroPrim_min )
            {
                data2 = data_vec2[i];
                is_found2 = true;
                found_index2=i;
            }

            double alpha = root.m_alpha;
            double radius = root.m_radius;
            if(is_found2)
            {
                if(data2.get_nalpha()<m_NAlpha || data2.get_nradius()<m_NRadius)
                {
                    data_vec2.at(found_index2)=DataS2(DeltaOmZeroPrim, ro, alpha, radius, m_NAlpha, m_NRadius);
                    dataS2_save(data_vec2, path);
                }
            }
            else{
                data_vec2.push_back(DataS2(DeltaOmZeroPrim, ro, alpha, radius, m_NAlpha, m_NRadius));
                sort(data_vec2.begin(), data_vec2.end(), compare_DataS2_omega0p);
                dataS2_save(data_vec2, path);
            }
        }
    }

    return base_prim_1vector;
}


std::vector<double> GridS2::grid(double c0, double xi, unsigned int nfft,
        bool s2use, bool s2save, std::string path) const
{

    std::vector<double> base_prim_1vector=grid_prim(c0, nfft, s2use, s2save, path);

    //double ro = density(base_prim_1vector);
    std::vector<double> basis(16, 0.0);
    double r_scale=sqrt(1.0-c0);

    //basis[0]=ro;
    std::vector<double> temp;
    temp = num::multiply_AB(base_prim_1vector, num::inverse( num::cholesky( m_fm->postrmf(xi), 4), 4), 4, 4, 4);
    for(int i=0; i<16; i++)
        basis[i]=r_scale*temp[i]; //basis[i+1]=

    return basis;
}

// convert vector from hyper-sphere space to hyper-ellipsoid space
std::vector<double> GridS2::convert(double c0, double xi, const std::vector<double>& grid_prim) const
{
    std::vector<double> ws1(16, 0.0);
    double r_scale=sqrt(1.0-c0);

    std::vector<double> temp(16, 0.0);
    temp = num::multiply_AB(grid_prim, num::inverse( num::cholesky( m_fm->postrmf(xi), 4), 4), 4, 4, 4);
    for(int i=0; i<16; i++)
        ws1[i]=r_scale*temp[i];

    return ws1;
}

std::vector<DataS2> GridS2::dataS2_load(std::string path) const
{
    double omega0p, density, alpha, radius;
    int n_alpha, n_radius;

	std::ifstream in(path.c_str());
	if(!in.is_open())
	{
        std::string error = "Can't open file: "+path;
        throw std::runtime_error(error);
	}

    std::vector<DataS2> t_data;
    t_data.reserve(100);
    std::string line;
    in >> std::skipws;
    std::getline(in,line);
    while(in >> std::skipws >> omega0p >> density >> alpha >> radius >> n_alpha >> n_radius){
        t_data.push_back(DataS2(omega0p, density, alpha, radius, n_alpha, n_radius));
    }
    in.close();

    return t_data;
}

void GridS2::dataS2_save(std::vector<DataS2>& data, std::string path="dataS2.txt") const
{
    std::string s1="DelOmZeroP";        // \Delta\omega_{0}^{'}
    std::string s2="Density";           // density (ro)
    std::string s3="Alpha";             // alpha (angle)
    std::string s4="Radius";             // covering radius
    std::string s5="NAlpha";            // incrementation deep in root finding algorithm 'find_alpha'
    std::string s6="NRadius";           // incrementation deep in algorithm obtaining radius of covering
    std::string tab="\t";
    std::string tab2="\t\t";
    std::string first_line=s1+tab+s2+tab2+s3+tab2+s4+tab2+s5+tab+s6+"\n";

    std::fstream fs;
    fs.open( path, std::fstream::out | std::fstream::trunc);
    if( fs.good() == true )
    {
        fs << std::left << first_line;// << std::fixed;
        for(unsigned int i = 0; i<data.size(); i++)     // std::fixed
        {
            fs  << std::setw(15) << std::setprecision(14) << data[i].m_omega0p << tab
                << std::setprecision(5) << data[i].m_density << tab2
                << std::setw(15) << std::setprecision(14) << data[i].m_alpha << tab
                << std::setw(14)<< std::setprecision(13) << data[i].m_radius << tab
                << std::setw(2) << data[i].get_nalpha() << tab
                << data[i].get_nradius() << "\n";
        }
        //fs.unget();
        fs.close();

    }
    else
    {
        std::string error = "Can't open: " + path;
        throw std::runtime_error(error);
    }
}

int GridS2::check_nalpha(int n_alpha)
{
    if(n_alpha<=0)
        n_alpha=1;

    return n_alpha;
}

int GridS2::check_nradius(int n_radius)
{
    if(n_radius<=0)
        n_radius=1;

    return n_radius;
}
