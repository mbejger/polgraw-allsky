///////////////////////////////////////
// Name:        Manual.h
// Author:      Andrzej Pisarski
// Copyright:   Andrzej Pisarski
// License:     CC-BY-NC-ND
// Created:     15/11/2015
// Modification:14/06/2020 A.Pisarski
///////////////////////////////////////

#ifndef MANUAL_H
#define MANUAL_H

#include <vector>
#include <string>

class Manual
{
    public:
        Manual(std::string, unsigned int pages=15);

        std::vector<std::string> help();
        std::vector<std::string> version();
        void save_help();
        void print_help();
        void print_version();
        void print_author();

    private:
        std::vector<std::string> m_help;
        std::vector<std::string> m_version;
        unsigned int m_pagesNr;    /// Number of pages

        std::vector<std::string> init_m_help() const;
        std::vector<std::string> init_m_version(std::string) const;
};

#endif // MANUAL_H
