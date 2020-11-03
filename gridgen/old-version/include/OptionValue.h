///////////////////////////////////////
// Name:        OptionValue.h
// Author:      Andrzej Pisarski
// Copyright:   Andrzej Pisarski
// License:     CC-BY-NC-ND
// Created:     20/05/2015
///////////////////////////////////////

#ifndef OPTIONVALUE_H
#define OPTIONVALUE_H

#include <string>

class OptionValue
{
    public:
        OptionValue();
        OptionValue(int, int, double, std::string);

        int m_k;
        int m_i;
        double m_d;
        std::string m_s;
};

#endif // OPTIONVALUE_H

/// OptionValue's members can not be 'const' - they should stay changeable (eg. in ReadConfig)
/// ?? private + friend class ReadConfig + get_options() ??

/// use const for reference and pointer parameters
/// It's about self-documenting your code and your assumptions.
/// If your code has many people working on it and your functions are non-trivial
/// then you should mark "const" any and everything that you can.
/// When writing industrial-strength code, you should always assume
/// that your coworkers are psychopaths
///
/// Dropping const from function prototypes has the advantage that you don't need to
/// alter the header file if you decide to drop const from implementation part later
/// description source:
/// http://stackoverflow.com/questions/117293/use-of-const-for-function-parameters
/// const - stala wartosc przez cale trwanie programu
/// static - dostep do pola nawet bez tworzenia obiektu danej klasy
/// References are inherently const, that is you can't change what they refer to.
/// references are always const.
