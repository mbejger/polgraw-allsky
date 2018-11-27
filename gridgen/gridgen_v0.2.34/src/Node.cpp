///////////////////////////////////////
// Name:        Node.cpp
// Author:      Andrzej Pisarski
// Copyright:   Andrzej Pisarski
// License:     CC-BY-NC-ND
// Created:     20/05/2015
///////////////////////////////////////

#include "Node.h"

Node::Node(): m_distance(0.), m_coef(0) {}
Node::Node(double distance, const std::vector<int>& coeff): m_distance(distance), m_coef(coeff) //<int>&
{
    //ctor
}

Node::Node(const Node& n): m_distance(n.m_distance), m_coef(n.get_coef())
{

}

void Node::operator=(const Node& n)
{
    m_distance=n.m_distance;
    m_coef=n.get_coef();
}

std::vector<int> Node::get_coef() const
{
    return m_coef;
}

bool compare_Node_distance(const Node& na, const Node& nb)
{
    return na.m_distance < nb.m_distance; //(na.m_distance - nb.m_distance) < 0.0;
}

bool compare_Node_distance_decrease(const Node& na, const Node& nb)
{
    return na.m_distance > nb.m_distance; //(na.m_distance - nb.m_distance) < 0.0;
}


/*bool Node::operator < (const Node& nb) const
{
    return m_distance < nb.m_distance;
}*/

