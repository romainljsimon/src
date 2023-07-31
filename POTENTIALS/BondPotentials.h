
#ifndef BONDPOTENTIALS_H_
#define BONDPOTENTIALS_H_

#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include "../INPUT/Parameter.h"

class BondPotentials
{

private:
    const int m_particleTypes {};
    const std::vector<double> m_bondPotentials {};

public:
    // Bonds constructor

    BondPotentials() = default;

    explicit BondPotentials (param::Parameter param)
    : m_particleTypes (param.get_int("particleTypes"))
    , m_bondPotentials(initializeBondPotentials(m_particleTypes))
    {}


    static std::vector<double> initializeBondPotentials(int particleTypes)
    {
        const int lenBonds { (particleTypes * (particleTypes + 1)) / 2 * 6};
        std::vector<double> bondPotentials(lenBonds);

        param::Parameter potentials("./potentials.txt");
        std::string keyBond{"bondCoeff"};
        std::string delimiter{"|"};

        for (int i = 1; i <= particleTypes; i++)
        {
            const std::string strI{std::to_string(i)};

            for (int j = i; j <= particleTypes; j++)
            {

                const std::string strJ{std::to_string(j)};
                std::string keyIJ{keyBond};
                keyIJ.append(strI).append(strJ);
                std::string bondCoeff{ potentials.get_string(keyIJ, "0.|0.|0.|0.|0.|0.|")};
                const int indexIJ {(j - i + particleTypes * (i - 1) - ((i - 2) * (i-1)) / 2) * 6};
                size_t pos;
                std::string token;
                std::vector<std::string> coeffList{};

                while ((pos = bondCoeff.find(delimiter)) != std::string::npos)
                {
                    token = bondCoeff.substr(0, pos);
                    coeffList.push_back(token);
                    bondCoeff.erase(0, pos + delimiter.length());
                }

                const double& r0{std::stod(coeffList[1])};
                bondPotentials[indexIJ + 0] = r0 * r0;// k constant for bond i


                bondPotentials[indexIJ + 1] = std::stod(coeffList[0]); // squareR0 constant for bond i

                bondPotentials[indexIJ + 3] = 4 * std::stod(coeffList[2]); // Epsilon constant for bond i

                const double& sigma{std::stod(coeffList[3])};
                bondPotentials[indexIJ + 4] = sigma * sigma; // square Sigma constant for bond i

                const double& rc{std::stod(coeffList[4])};

                bondPotentials[indexIJ + 2] = rc * rc; // Rc constant for bond i

                bondPotentials[indexIJ + 5] = std::stod(coeffList[5]); // Shift constant for bond i

                //std::cout << i << " "<< j << " "<< indexIJ << " "<< bondPotentials[indexIJ+0] << " "
                //<< bondPotentials[indexIJ+1] << " "<< bondPotentials[indexIJ+2] << " "<< bondPotentials[indexIJ+3] << "\n" ;

            }
        }

        return bondPotentials;
    }


    [[nodiscard]] double feneBondEnergyIJ(const double &squareDistance, const int &particleTypeI,
                                          const int &particleTypeJ) const;

    //[[nodiscard]] std::vector<double> getPotentialsIJ(const int &i, const int &j) const;

    [[nodiscard]] int getIndexIJ(const int &i, const int &j) const;

};
#endif /* BONDPOTENTIALS_H_ */
