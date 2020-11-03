///////////////////////////////////////
// Name:        ReadConfig.h
// Author:      Andrzej Pisarski
// Copyright:   Andrzej Pisarski
// License:     CC-BY-NC-ND
// Created:     20/05/2015
///////////////////////////////////////

#ifndef READCONFIG_H
#define READCONFIG_H

#include "OptionValue.h"
#include <vector>
#include <string>
//C11
class ReadConfig
{
    public:
        ReadConfig();
        ReadConfig(const std::string&);
        ReadConfig(const std::vector<OptionValue>);///?
        ~ReadConfig();// delete?

        std::vector<std::string> get_options() const;
        std::vector<OptionValue> get_values() const;

        //int flag;           // file stream state

    private:
        static const std::vector<std::string> m_options;    /// const? // vector of available options
        std::vector<OptionValue> m_values;                  // storage values of m_options (order is important!)

        bool set_line(std::string::const_iterator, std::string::const_iterator) const;
        std::vector<OptionValue> set_command(const std::vector<std::string>&) const;
        std::vector<OptionValue> read_config_file(const std::string& path="gg.ini") const;   // read option's values from file 'gg.ini'
};

#endif // READCONFIG_H
