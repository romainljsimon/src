#ifndef MOLECULES_H_
#define MOLECULES_H_

#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <vector>
#include <cmath>
#include <algorithm>
#include <typeinfo>
#include <tuple>
#include "../INPUT/Parameter.h"
#include "../POTENTIALS/PairPotentials.h"
#include "../POTENTIALS/BondPotentials.h"


class Neighbors;

class Molecules
{
    friend class Neighbors;

private:
    // POTENTIALS constructor
    const int m_nDims {3};
    const PairPotentials m_systemPairPotentials {};
    const BondPotentials m_systemBondPotentials {};
    const int m_nParticles {};
    const double m_lengthCube {};
    const double m_halfLengthCube {};
    const std::vector<int> m_bondsArray {};
    const std::vector<int> m_bondsIndex {};
    std::vector<int> m_newFlags;
    std::vector<int> m_flagsArray {};
    std::vector<double> m_positionArray {};
    std::vector<int> m_particleTypeArray {};
    std::vector<int> m_moleculeTypeArray {};
    const std::string m_saveHeaderString{};
    using PosIterator = std::vector<double>::const_iterator;
    using BondsIterator = std::vector<int>::const_iterator;

    //std::vector<double> typeArray{1, 2};
    Molecules (param::Parameter param, PairPotentials  systemPairPotentials,
               BondPotentials  systemBondPotentials, const std::string& path,
               const std::tuple<std::vector<int>, std::vector<int>>& bondsData)

            : m_nParticles ( initializeNumber(path))
            , m_systemPairPotentials(std::move(systemPairPotentials))
            , m_systemBondPotentials(std::move(systemBondPotentials))
            , m_lengthCube {pow (static_cast<double>(m_nParticles) / param.get_double("density"), 1. / 3.) }
            , m_halfLengthCube (0.5 * m_lengthCube)
            , m_bondsArray(std::get<0>(bondsData))
            , m_bondsIndex(std::get<1>(bondsData))
            , m_saveHeaderString(initializeHeaderString(m_nParticles, m_lengthCube))

    {
        m_flagsArray.resize(m_nDims * m_nParticles, 0);
        initializeParticles(path);

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

        int col { 8 };

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

                else if (c < 5)
                {
                    infile >> m_positionArray[m_nDims * r + c - (col - 6)];
                    //Take INPUT from file and put into positionArray
                }
                else
                {
                    infile >> m_flagsArray[m_nDims * r + c - (col - 3)];
                }
            }
        }
        for (auto it=m_positionArray.begin(); it < m_positionArray.end(); it+=m_nDims)
        {
            periodicBC(it);
            reinitializeFlags();
        }
        infile.close();

    }

public:
    Molecules (param::Parameter param, PairPotentials  systemPairPotentials,
               BondPotentials  systemBondPotentials, const std::string& path)

        : Molecules(param, std::move(systemPairPotentials), std::move(systemBondPotentials), path, initializeBondsArray(initializeNumber(path)))
    {
    }

    static std::string initializeHeaderString(const int& nParticles, const double& lengthCube)
    {
        std::string headerString {std::to_string(nParticles)};
        headerString.append("\n");
        std::string lengthStr {std::to_string(lengthCube)};
        std::string zeroString {" 0.0 0.0 0.0 "};
        headerString.append("Lattice=");
        headerString.append("\"");
        headerString.append(lengthStr);
        headerString.append(zeroString);
        headerString.append(lengthStr);
        headerString.append(zeroString);
        headerString.append(lengthStr);
        headerString.append("\"");
        headerString.append(" Properties=molecule_type:S:1:type:I:1:pos:R:3:");
        headerString.append("\n");
        return headerString;
    }

    static std::tuple<std::vector<int>, std::vector<int>> initializeBondsArray(const int& nParticles)
    {
        std::ifstream infile("./bonds.txt");

        if (!infile.is_open())
            std::cout << "Error opening file";

        int nBonds{};
        infile >> nBonds;


        std::vector<int> bondsArray(2 * nBonds);
        std::vector<int> bondsIndex(nParticles + 1, 0);

        for (int r = 0; r < nBonds; r++) //Outer loop for rows
        {

            int indexI {};
            infile >> indexI;

            int indexJ {};
            infile >> indexJ;

            std::rotate(bondsArray.rbegin(),
                        bondsArray.rbegin() + 1, bondsArray.rend() - bondsIndex[indexI+1] );
            bondsArray[bondsIndex[indexI+1]] = indexJ;
            std::for_each(bondsIndex.begin() + indexI+1, bondsIndex.end(), [](int& d) { d+=1;});

            std::rotate(bondsArray.rbegin(), bondsArray.rbegin(), bondsArray.rend()-bondsIndex[indexI+1] );
            bondsArray[bondsIndex[indexJ+1]] = indexI;
            std::for_each(bondsIndex.begin() + indexJ+1, bondsIndex.end(), [](int& d) { d+=1;});
            //bondsArray.insert(bondsArray.begin() + bondsIndex[indexI], indexJ);
            //bondsArray[indexI].push_back(indexJ);
            //bondsArray[indexJ].push_back(indexI);

        }
        infile.close();
        /***

        std::cout << bondsArray.size() << "\n";
         ***/
        return std::make_tuple(bondsArray, bondsIndex);
    }

    void updateFlags(const int& indexParticle);

    template<typename InputIt>
    void updatePositionI(const int& i, InputIt newPosItBegin)
    {
        updatePositionI(i, newPosItBegin, m_nDims);
    }

    template<typename InputIt>
    void updatePositionI(const int& i, InputIt newPosItBegin, const int& lenPos)
    {
        const int realIndex {i * m_nDims};
        std::copy(newPosItBegin, newPosItBegin + lenPos, m_positionArray.begin() + realIndex);
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

        const auto &bondsItBegin { getBondsItBeginI(indexParticle) };
        const auto &bondsItEnd {getBondsItEndI(indexParticle)};
        const int& particleTypeI { m_particleTypeArray[indexParticle] };

        for (auto it = bondsItBegin; it < bondsItEnd; it++)
        {
            const int& indexJ {*it};
            const int& particleTypeJ { m_particleTypeArray[indexJ] };

            const double squareDistance { squareDistancePair(posItBegin,
                                                                 getPosItBeginI(indexJ)) };
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

        for (auto it = NeighItBegin; it < NeighItBegin + lenNeigh; ++it)
        {
            const int& indexJ {*it};

            const int& typeJ { m_particleTypeArray[indexJ] };
            const double squareDistance { squareDistancePair(posItBegin,
                                                             getPosItBeginI(indexJ)) };
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

        const auto& posItBegin {getPosItBeginI(indexParticle)};
        return energyParticleMolecule(indexParticle, posItBegin, NeighItBegin, lenNeigh);
    }

    void reinitializeFlags();

    template<typename InputPosIt, typename InputNeighIt>
    double energyPairParticleExtraMolecule(const int& indexParticle, InputPosIt posItBegin,
                                  InputNeighIt NeighItBegin, const int& lenNeigh, const int& typeMoleculeI) const
    {

        double energy { 0. };
        const int& particleType {m_particleTypeArray[indexParticle]};

        for (auto it = NeighItBegin; it < NeighItBegin + lenNeigh; ++it)
        {
            const int indexJ {*it};
            const int& typeMoleculeJ {m_moleculeTypeArray[indexJ]};

            if (typeMoleculeJ != typeMoleculeI)
            {
                const int& typeJ {m_particleTypeArray[indexJ]};
                const double squareDistance { squareDistancePair(posItBegin, m_positionArray.begin() + m_nDims * indexJ)};
                energy += m_systemPairPotentials.ljPairEnergy(squareDistance, particleType, typeJ);

            }
        }
        return energy;

    }

    template<typename InputNeighIt>
    double energyPairParticleExtraMolecule(const int& indexParticle,
                                       InputNeighIt NeighItBegin, const int& lenNeigh, const int& typeMoleculeI) const
    {
        auto posItBegin {m_positionArray.begin() + m_nDims * indexParticle};
        return energyPairParticleExtraMolecule(indexParticle, posItBegin, NeighItBegin, lenNeigh, typeMoleculeI);
    }



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

        const auto &bondsItBegin { getBondsItBeginI(indexParticle) };
        const auto &bondsItEnd {getBondsItEndI(indexParticle)};
        const int& particleTypeI { m_particleTypeArray[indexParticle] };
        const int& swapParticleTypeI {m_particleTypeArray[indexSwap]};

        for (auto it = bondsItBegin; it < bondsItEnd; it++)
        {
            const int& indexJ {*it};

            if (indexJ != indexSwap)
            {
                const int& particleTypeJ { m_particleTypeArray[indexJ] };

                const double squareDistance { squareDistancePair(posItBegin,
                                                                 m_positionArray.begin() + m_nDims*indexJ) };

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

    [[nodiscard]] PosIterator getPosItEndI(const int &i) const
    {
        const int realIndex { (i+1) * m_nDims };
        auto outIt {m_positionArray.begin() + realIndex};
        return outIt;
    }

    [[nodiscard]] BondsIterator getBondsItBeginI(const int &i) const
    {
        auto outIt {m_bondsArray.begin() + m_bondsIndex[i]};
        return outIt;
    }

    [[nodiscard]] BondsIterator getBondsItEndI(const int &i) const
    {
        auto outIt {m_bondsArray.begin() + m_bondsIndex[i + 1]};
        return outIt;
    }


    [[nodiscard]] const int& getNDims() const;

    [[nodiscard]] double getCosAngleMolecule(const std::vector<int>& orderVector, const int& indexMolecule) const;

    [[nodiscard]] std::vector<int> getOrderVector(const int &indexMolecule) const;

    int getNMolecules() const;
};

#endif /* MOLECULES_H_ */
