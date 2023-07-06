/*
MIT License
Copyright (c) 2019 - present H. Watanabe
The latest version is available at
https://github.com/kaityo256/params

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NON INFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#ifndef BONDS_H_
#define BONDS_H_

#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include "Parameter.h"

class Bonds
{

private:
    std::vector<double> m_bondKIArray {};
    std::vector<double> m_bondR0IArray {};
    std::vector<double> m_bondEpsilonIArray {};
    std::vector<double> m_bondSigmaIArray {};
    std::vector<double> m_bondRcIArray {};
    const int m_bondTypes {};


public:
    // Polymer constructor
    explicit Bonds (param::Parameter param)
    : m_bondTypes (param.get_int( "bondTypes"))
    {
        // Bonds Initializations
        initializeBonds();
    }


    void initializeBonds()
    {
        param::Parameter potentials("./potentials.txt" );
        std::string keyBond {"bondCoeff"};
        std::string delimiter { "|"};
        for (int i=0; i<m_bondTypes; i++)
        {
            std::string strI { std::to_string(i+1)};
            std::string keyI {keyBond};
            keyI.append(strI);
            std::string potentialCoeff {potentials.get_string(keyI)};
            size_t pos = 0;
            std::string token;
            std::vector<std::string> coeffList {};

            while ((pos = potentialCoeff.find(delimiter)) != std::string::npos)
            {
                token = potentialCoeff.substr(0, pos);
                coeffList.push_back(token);
                potentialCoeff.erase(0, pos + delimiter.length());
            }
            m_bondKIArray[i] = std::stod(coeffList[0]);
            m_bondR0IArray[i] = std::stod(coeffList[1]);
            m_bondEpsilonIArray[i] = std::stod(coeffList[2]);
            m_bondSigmaIArray[i] = std::stod(coeffList[3]);
            m_bondRcIArray[i] = std::stod(coeffList[4]);

        }
    }
};
#endif /* BONDS_H_ */
