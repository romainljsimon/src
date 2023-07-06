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

#ifndef POTENTIALS_H_
#define POTENTIALS_H_

#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include "Parameter.h"

class Potentials
{

private:

    std::vector<double> m_epsilonIJArray {};
    std::vector<double> m_sigmaIJArray {};
    std::vector<double> m_rcIJArray {};
    const int m_particleTypes {};




public:
    // Polymer constructor
    explicit Potentials (param::Parameter param)
    : m_particleTypes (param.get_int( "particleTypes"))
    {
        // Potential Initializations
        initializePotentials();
    }


    void initializePotentials()
    {

        int lenPairs {m_particleTypes * (m_particleTypes + 1) / 2};
        m_epsilonIJArray.resize(lenPairs);
        m_sigmaIJArray.resize(lenPairs);
        m_rcIJArray.resize(lenPairs);
        param::Parameter potentials("./potentials.txt" );
        std::string keyPair {"pairCoeff"};
        std::string delimiter { "|"};
        for (int i=1; i<=m_particleTypes; i++)
        {
            std::string strI { std::to_string(i)};
            for (int j=i; j<=m_particleTypes; j++)
            {
                std::string strJ { std::to_string(j)};
                std::string keyIJ {keyPair};
                keyIJ.append(strI).append(strJ);
                std::string potentialCoeff {potentials.get_string(keyIJ)};
                const int indexIJ {j - i + m_particleTypes * (i - 1) - ((i - 1) * i) / 2};
                size_t pos = 0;
                std::string token;
                std::vector<std::string> coeffList {};

                while ((pos = potentialCoeff.find(delimiter)) != std::string::npos) {
                    token = potentialCoeff.substr(0, pos);
                    coeffList.push_back(token);
                    potentialCoeff.erase(0, pos + delimiter.length());
                }
                m_epsilonIJArray[indexIJ] = std::stod(coeffList[0]);
                m_sigmaIJArray[indexIJ] = std::stod(coeffList[1]);
                m_rcIJArray[indexIJ] = std::stod(coeffList[2]);

            }
        }
    }
};
 // namespace param
#endif /* POTENTIALS_H_ */
