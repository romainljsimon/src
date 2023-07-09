
#ifndef BONDPOTENTIALS_H_
#define BONDPOTENTIALS_H_

#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include "../INPUT/Parameter.h"
#include "../PARTICLES/Particles.h"

class BondPotentials
{

private:
    const int m_particleTypes {};
    const std::vector<std::vector<double>> m_bondPotentials {};
    const std::vector<std::vector<int>> m_bondsArray {};


public:
    // Bonds constructor

    BondPotentials() = default;

    explicit BondPotentials (param::Parameter param, const std::string& path)
    : m_particleTypes (param.get_int("particleTypes"))
    , m_bondPotentials(initializeBondPotentials(m_particleTypes))
    , m_bondsArray(initializeBondsArray(path))
    {}


    static std::vector<std::vector<double>> initializeBondPotentials(int particleTypes) {

        std::vector<std::vector<double>> bondPotentials((particleTypes * (particleTypes - 1)) / 2,
                                                        std::vector<double>(6));
        param::Parameter potentials("./POTENTIALS.txt");
        std::string keyBond{"bondCoeff"};
        std::string delimiter{"|"};

        for (int i = 1; i < particleTypes; i++)
        {
            std::string strI{std::to_string(i)};

            for (int j = i + 1; j <= particleTypes; j++)
            {

                std::string strJ{std::to_string(j)};
                std::string keyIJ{keyBond};
                keyIJ.append(strI).append(strJ);
                std::string bondCoeff{potentials.get_string(keyIJ, "0.|0.|0.|0.|0.|0.|")};
                const int indexIJ{j - i - 1 + particleTypes * (i - 1) - ((i - 1) * i) / 2};
                size_t pos = 0;
                std::string token;
                std::vector<std::string> coeffList{};

                while ((pos = bondCoeff.find(delimiter)) != std::string::npos)
                {
                    token = bondCoeff.substr(0, pos);
                    coeffList.push_back(token);
                    bondCoeff.erase(0, pos + delimiter.length());
                }

                bondPotentials[indexIJ][0] = std::stod(coeffList[0]);// k constant for bond i

                double r0{std::stod(coeffList[1])};
                bondPotentials[indexIJ][1] = r0 * r0; // squareR0 constant for bond i

                bondPotentials[indexIJ][2] = std::stod(coeffList[2]); // Epsilon constant for bond i

                double sigma{std::stod(coeffList[3])};
                bondPotentials[indexIJ][3] = sigma * sigma; // square Sigma constant for bond i

                double rc{std::stod(coeffList[4])};

                bondPotentials[indexIJ][4] = rc * rc; // Rc constant for bond i

                bondPotentials[indexIJ][5] = std::stod(coeffList[5]); // Shift constant for bond i
            }
        }

        return bondPotentials;
    }

    static std::vector<std::vector<int>> initializeBondsArray(const std::string& path)
    {
        std::ifstream infile(path);

        if (!infile.is_open())
            std::cout << "Error opening file";


        int nParticles{};
        infile >> nParticles;

        int nBonds{};
        infile >> nBonds;

        std::vector<std::vector<int>> bondsArray(nParticles);

        for (int r = 0; r < nBonds; r++) //Outer loop for rows
        {

            int indexI {};
            infile >> indexI;

            int indexJ {};
            infile >> indexJ;

            bondsArray[indexI].push_back(indexJ);
            bondsArray[indexJ].push_back(indexI);

        }
        infile.close();
        /***

        std::cout << bondsArray.size() << "\n";
         ***/
        return bondsArray;
    }




    [[nodiscard]] double feneBondEnergyIJ(const double &squareDistance, const int &particleTypeI,
                                          const int &particleTypeJ) const;

    [[nodiscard]] std::vector<double> getPotentialsIJ(const int &i, const int &j) const;

    [[nodiscard]] double feneBondEnergyI(const int &indexParticle, const std::vector<double> &positionParticle,
                                         const Particles& systemParticles, const int &indexSkip) const;

    [[nodiscard]] int getIndexIJ(const int &i, const int &j) const;

    [[nodiscard]] std::vector<int> getBondsI(const int &i) const;
};
#endif /* BONDPOTENTIALS_H_ */
