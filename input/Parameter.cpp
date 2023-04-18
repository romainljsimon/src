/*
 * param.cpp
 *
 *  Created on: 17 apr. 2023
 *      Author: Romain Simon
 */
#include <string>
#include "Parameter.h"



template <>
bool param::Parameter::get(std::string key)
{
    check_key(key);
    return get<bool>(key, false);
}

template <>
int param::Parameter::get(std::string key, int value)
{
    if (!contains(key))
    {
        return value;
    }
    return std::stoi(params[key]);
}

template <>
int param::Parameter::get(std::string key)
{
    check_key(key);
    return get<int>(key, 0);
}

template <>
double param::Parameter::get(std::string key, double value)
{
    if (!contains(key)) {
        return value;
    }
    return std::stod(params[key]);
}

template <>
double param::Parameter::get(std::string key)
{
    check_key(key);
    return get<double>(key, 0.0);
}
