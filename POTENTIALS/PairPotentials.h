#ifndef POTENTIALS_H_
#define POTENTIALS_H_

#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include "INPUT/Parameter.h"
#include "PARTICLES/Particles.h"

class PairPotentials
{

private:
    const int m_nParticleTypes {};
    const std::vector<double> m_pairPotentials {};



public:
    // POTENTIALS constructor

    explicit PairPotentials (param::Parameter param)
    : m_nParticleTypes (param.get_int( "particleTypes"))
    , m_pairPotentials (initializePotentials(m_nParticleTypes))
    {
    }


    static std::vector<double> initializePotentials(int nParticleTypes)
    {

        int lenPairs {nParticleTypes * (nParticleTypes + 1) / 2};
        std::vector<double> pairPotentials (lenPairs * 4) ; //, std::vector<double>(4));

        param::Parameter potentials("./POTENTIALS.txt" );
        std::string keyPair {"pairCoeff"};
        std::string delimiter { "|"};
        for (int i=1; i<=nParticleTypes; i++)
        {
            std::string strI { std::to_string(i)};

            for (int j=i; j<=nParticleTypes; j++)
            {
                std::string strJ { std::to_string(j)};
                std::string keyIJ {keyPair};
                keyIJ.append(strI).append(strJ);
                std::string potentialCoeff {potentials.get_string(keyIJ)};
                const int indexIJ {j - i + nParticleTypes * (i - 1) - ((i - 2) * (i-1)) / 2};
                size_t pos;
                std::string token;
                std::vector<std::string> coeffList {};

                while ((pos = potentialCoeff.find(delimiter)) != std::string::npos) {
                    token = potentialCoeff.substr(0, pos);
                    coeffList.push_back(token);
                    potentialCoeff.erase(0, pos + delimiter.length());
                }


                pairPotentials[indexIJ * 4 + 0] = std::stod(coeffList[0]); // Epsilon IJ pair constant
                const double& sigmaIJ {std::stod(coeffList[1])};
                pairPotentials[indexIJ * 4 + 1] = sigmaIJ * sigmaIJ; // SigmaSquare IJ pair constant
                const double& rcIJ { std::stod(coeffList[2]) };
                pairPotentials[indexIJ * 4 + 2] = rcIJ * rcIJ; // RcSquare IJ pair constant
                pairPotentials[indexIJ * 4 + 3] = std::stod(coeffList[3]); // Shift IJ
                /***
                std::cout << pairPotentials[indexIJ * 4 + 0] << " ";
                std::cout << pairPotentials[indexIJ * 4 + 1] << " ";
                std::cout << pairPotentials[indexIJ * 4 + 2] << " ";
                std::cout << pairPotentials[indexIJ * 4 + 3] << " ";
                std::cout << "\n";
                ***/

            }
        }

        return pairPotentials;
    }

    std::vector<double> getPotentialsIJ(const int& i, const int& j, const int& k) const;

    [[nodiscard]] std::vector<double> getPotentialsIJ(const int &i, const int &j) const;

    [[nodiscard]] double ljPairEnergy(const double &squareDistance, const int &typeI, const int &typeJ) const;

    [[nodiscard]] double energyPairParticle(const int &indexParticle, const std::vector<double> &positionParticle,
                                            const Particles &systemParticles, const std::vector<int> &neighborIList,
                                            const int &indexSkip) const;

    [[nodiscard]] int getParticleTypes() const;

    [[nodiscard]] int getIndexIJ(const int &i, const int &j) const;

    [[nodiscard]] double getSquareRcIJ(const int &i, const int &j) const;
};
 // namespace param
#endif /* POTENTIALS_H_ */
