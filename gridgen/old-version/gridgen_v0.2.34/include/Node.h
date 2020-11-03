///////////////////////////////////////
// Name:        Node.h
// Author:      Andrzej Pisarski
// Copyright:   Andrzej Pisarski
// License:     CC-BY-NC-ND
// Created:     20/05/2015
///////////////////////////////////////

// Node class storage two values:
// distance from origin of coordinate axis (double; hyperspace with hyper-sphere)
// and node given by vector of 4 coefficients (vector<int>; (scalars)),
// (4 coefficients and base vectors crate linear combination).

#ifndef Node_H
#define Node_H

#include<vector>

class Node
{
    public:
        Node();
        Node(double, const std::vector<int>&);///&
        Node(const Node&);

        void operator=(const Node&);

        std::vector<int> get_coef() const;
        //bool operator < (const Node& ) const;
        double m_distance;

    private:
        std::vector<int> m_coef;
};

bool compare_Node_distance(const Node&, const Node&);
bool compare_Node_distance_decrease(const Node&, const Node&);

#endif // Node_H
