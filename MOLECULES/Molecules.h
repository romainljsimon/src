#ifndef MOLECULES_H_
#define MOLECULES_H_

#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <vector>
#include <cmath>
#include <cmath>
#include <typeinfo>
#include "../INPUT/Parameter.h"
#include "../POTENTIALS/PairPotentials.h"
#include "../POTENTIALS/BondPotentials.h"


class Neighbors;

class Molecules
{
    friend class Neighbors;

public:
    // POTENTIALS constructor
    const int m_nDims {3};
    const PairPotentials m_systemPairPotentials {};
    const BondPotentials m_systemBondPotentials {};
    const std::vector<std::vector<int>> m_bondsArray {};
    const int m_nParticles {};
    const double m_lengthCube {};
    const double m_halfLengthCube {};
    std::vector<int> m_newFlags;
    std::vector<int> m_flagsArray {};
    std::vector<double> m_positionArray {};
    std::vector<int> m_particleTypeArray {};
    std::vector<int> m_moleculeTypeArray {};
    const std::string m_saveHeaderString{};
    using PosIterator = std::vector<double>::const_iterator;

    //std::vector<double> typeArray{1, 2};



    Molecules (param::Parameter param, PairPotentials  systemPairPotentials,
               BondPotentials  systemBondPotentials, const std::string& path)

        : m_nParticles ( initializeNumber(path))
        , m_systemPairPotentials(std::move(systemPairPotentials))
        , m_systemBondPotentials(std::move(systemBondPotentials))
        , m_lengthCube {pow (static_cast<double>(m_nParticles) / param.get_double("density"), 1. / 3.) }
        , m_halfLengthCube (0.5 * m_lengthCube)
        , m_bondsArray(initializeBondsArray())
        , m_saveHeaderString(initializeHeaderString(m_nParticles, m_lengthCube))
    {
        m_flagsArray.resize(3 * m_nParticles, 0);
        initializeParticles(path);
        std::cout << typeid(m_positionArray.begin()).name();

    }

    static std::string initializeHeaderString(const int& nParticles, const double& lengthCube)
    {
        std::string headerString {static_cast<char>(nParticles)};
        headerString.append("\n");
        std::string lengthStr {std::to_string(lengthCube)};
        std::string zeroString {" 0.0 0.0 0.0 "};
        headerString.append("Lattice=");
        headerString.append(std::to_string('"'));
        headerString.append(lengthStr);
        headerString.append(zeroString);
        headerString.append(lengthStr);
        headerString.append(zeroString);
        headerString.append(std::to_string('"'));
        headerString.append(" Properties=molecule_type:S:1:type:I:1:pos:R:3:");
        headerString.append("\n");
        return headerString;
    }

    static std::vector<std::vector<int>> initializeBondsArray()
    {
        std::ifstream infile("./bonds.txt");

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

    static int initializeNumber(const std::string& path)
    {
        std::ifstream infile(path);

        if (!infile.is_open())
            std::cout << "Error opening file" ;

        int nParticles {};
        infile >> nParticles;
        infile.close();
        return nParticles;
    }

    void initializeParticles(const std::string& path)
    {
        std::ifstream infile(path);

        if (!infile.is_open())
            std::cout << "Error opening file" ;

        int row{};
        infile >> row;

        int col { 5 };

        std::string line;
        getline(infile, line);
        getline(infile, line);

        m_moleculeTypeArray.resize(row);
        m_particleTypeArray.resize(row);
        m_positionArray.resize(m_nDims * row);
        //std::vector<std::vector <double>> positionArray(row, std::vector<double>(3));
        //std::vector<double> typeArray (row);
        //std::vector<int> moleculeType (row , 1);

        //Defining the loop for getting INPUT from the file

        for (int r = 0; r < row; r++) //Outer loop for rows
        {

            for (int c = 0; c < col; c++) //inner loop for columns
            {

                if (c == 0)
                {
                    infile >> m_moleculeTypeArray[r];
                }
                else if (c == 1)
                {

                    infile >> m_particleTypeArray[r];
                }

                else
                {
                    infile >> m_positionArray[m_nDims * r + c - (col - 3)];
                     //Take INPUT from file and put into positionArray
                }
            }
        }
        infile.close();
    }

    void updateFlags(const int& indexParticle);

    template<typename InputIt>
    void updatePositionI(const int& i, InputIt newPosItBegin)
    {
        const int realIndex {i * m_nDims};
        std::copy(newPosItBegin, newPosItBegin + m_nDims, m_positionArray.begin() + realIndex);
    }


    template<typename InputIt>
    void periodicBC(InputIt posItBegin)
/*
 *This function is an implementation of the periodic Boundary conditions.
 *If a particle gets out of the simulation box from one of the sides, it gets back in the box from the opposite side.
 */
    {
        for (int i = 0; i < m_nDims; i++)
        {
            //std::cout <<  m_newFlags[i] << "\n";
            const auto& posI {*posItBegin};
            if (posI < 0)
            {
                *posItBegin += m_lengthCube;
                m_newFlags.push_back(-1);
            }
            else if (posI > m_lengthCube)
            {
                *posItBegin -= m_lengthCube;
                m_newFlags.push_back(1);
            }
            else
            {
                m_newFlags.push_back(0);
            }
            posItBegin++;
        }
    }

    [[nodiscard]] const int& getNParticles() const;

    [[nodiscard]]  std::vector<double> getPositionI(const int& i) const;


    [[nodiscard]] const int& getParticleTypeI(const int &i) const;

    [[nodiscard]] const int& getMoleculeTypeI(const int &i) const;

    void swapParticleTypesIJ(const int &i, const int &j, const int &typeI, const int &typeJ);

    [[nodiscard]] const double& getLengthCube() const;

    [[nodiscard]] const double& getHalfLengthCube() const;

    void saveInXYZ(const std::string& path) const;

    void swapParticleTypesIJ(const int &i, const int &j);

    [[nodiscard]] double energySystemMolecule(const Neighbors &systemNeighbors) const;


    template<typename InputIt>
    double feneBondEnergyI(const int& indexParticle, InputIt posItBegin) const
    {
        double energy { 0. };

        const std::vector<int>& bondsI { m_bondsArray[indexParticle] };
        const int& particleTypeI { m_particleTypeArray[indexParticle] };

        for (int const &indexJ: bondsI)
        {

            const double squareDistance { squareDistancePair(posItBegin,
                                                                 m_positionArray.begin() + m_nDims*indexJ) };
            const int& particleTypeJ { m_particleTypeArray[indexJ] };
            energy += m_systemBondPotentials.feneBondEnergyIJ(squareDistance, particleTypeI, particleTypeJ);

        }
        return energy;
    }


    template<typename InputPosIt, typename InputNeighIt>
    double energyPairParticle(const int& indexParticle, InputPosIt posItBegin,
                              InputNeighIt NeighItBegin, const int& lenNeigh) const
    {
        double energy { 0. };
        const int& particleType {m_particleTypeArray[indexParticle]};

        for (auto it = NeighItBegin; it < NeighItBegin + lenNeigh; it++)
        {
            const int& indexJ {*it};

            const int& typeJ { m_particleTypeArray[indexJ] };
            const double squareDistance { squareDistancePair(posItBegin,
                                                             m_positionArray.begin() + m_nDims*indexJ) };
            energy += m_systemPairPotentials.ljPairEnergy(squareDistance, particleType, typeJ);


        }
        return energy;
    }

    template<typename InputPosIt, typename InputNeighIt>
    double energyParticleMolecule(const int& indexParticle, InputPosIt posItBegin,
                                  InputNeighIt NeighItBegin, const int& lenNeigh) const
    {
        double energy { 0. };
        energy += energyPairParticle(indexParticle, posItBegin, NeighItBegin, lenNeigh);
        energy += feneBondEnergyI(indexParticle, posItBegin);
        return energy;
    }

    template<typename InputNeighIt>
    double energyParticleMolecule(int indexParticle, InputNeighIt NeighItBegin, const int& lenNeigh) const
    {

        const auto& posItBegin {m_positionArray.begin() + m_nDims * indexParticle};
        return energyParticleMolecule(indexParticle, posItBegin, NeighItBegin, lenNeigh);
    }

    void reinitializeFlags();

    [[nodiscard]] double energyPairParticleExtraMolecule(const int& indexParticle, const std::vector<double>& positionParticle,
                                                      const std::vector<int>& neighborIList, const int& typeMoleculeI) const;

    [[nodiscard]] double energyPairParticleExtraMolecule(const int& indexParticle, const std::vector<int>& neighborIList, const int& typeMoleculeI) const;

    template<typename InputPosIt, typename InputNeighIt>
    double energyPairParticleSwap(const int& indexParticle, InputPosIt posItBegin,
                                  InputNeighIt NeighItBegin, const int& lenNeigh,
                                  const int& indexSwap) const
    {
        double energy { 0. };
        const int& particleType {m_particleTypeArray[indexParticle]};
        const int& swapParticleType {m_particleTypeArray[indexSwap]};

        for (auto it = NeighItBegin; it < NeighItBegin + lenNeigh; it++)
        {
            const int& indexJ {*it};
            if (indexJ != indexSwap)
            {

                const int& typeJ { m_particleTypeArray[indexJ] };
                const double squareDistance { squareDistancePair(posItBegin,
                                                                 m_positionArray.begin() + m_nDims*indexJ) };

                energy += m_systemPairPotentials.ljPairEnergy(squareDistance, swapParticleType, typeJ);
                energy -= m_systemPairPotentials.ljPairEnergy(squareDistance, particleType, typeJ);

            }
        }
        return energy;
    }

    template<typename InputIt>
    [[nodiscard]] double feneBondEnergyISwap(const int& indexParticle, InputIt posItBegin,
                               const int& indexSwap) const
    {
        double energy { 0. };

        const std::vector<int>& bondsI { m_bondsArray[indexParticle] };
        const int& particleTypeI { m_particleTypeArray[indexParticle] };
        const int& swapParticleTypeI {m_particleTypeArray[indexSwap]};

        for (int const &indexJ: bondsI)
        {

            if (indexJ != indexSwap)
            {
                const double squareDistance { squareDistancePair(posItBegin,
                                                                 m_positionArray.begin() + m_nDims*indexJ) };
                const int& particleTypeJ { m_particleTypeArray[indexJ] };
                energy += m_systemBondPotentials.feneBondEnergyIJ(squareDistance,
                                                                  swapParticleTypeI, particleTypeJ);
                energy -= m_systemBondPotentials.feneBondEnergyIJ(squareDistance, particleTypeI, particleTypeJ);
            }
        }
        return energy;
    }
    template<typename InputNeighIt>
    double energyParticleMoleculeSwap(const int& indexParticle, InputNeighIt NeighItBegin, const int& lenNeigh,
                                      const int& indexSwap) const
    {
        double energy { 0. };
        const auto& posItBegin {m_positionArray.begin() + 3 * indexParticle};
        energy += energyPairParticleSwap(indexParticle, posItBegin, NeighItBegin, lenNeigh, indexSwap);
        energy += feneBondEnergyISwap(indexParticle, posItBegin, indexSwap);
        return energy;
    }


    template<typename InputItI, typename InputItJ>
    double squareDistancePair(InputItI firstI, InputItJ firstJ) const
    {
        double squareDistance { 0. };

        for (int i = 0; i < m_nDims; i++)
        {
            double diff = *firstI - *firstJ;
            double absDiff {fabs(diff) };
            if (absDiff> m_halfLengthCube)
            {
                if (diff < 0.0)
                {
                    diff += m_lengthCube;
                }
                else
                {
                    diff -= m_lengthCube;
                }
            }

            squareDistance += diff * diff;
            firstI++;
            firstJ++;
        }

        return squareDistance;
    }

    [[nodiscard]] PosIterator getPosItBeginI(const int &i) const
    {
        const int realIndex { i * m_nDims };
        auto outIt {m_positionArray.begin() + realIndex};
        return outIt;
    }



};

#endif /* MOLECULES_H_ */
