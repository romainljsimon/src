/*
 * MonteCarlo.h
 *
 *  Created on: 27 oct. 2022
 *      Author: Romain Simon
 */

#ifndef NEIGHBORS_H_
#define NEIGHBORS_H_

#include <utility>
#include <iterator>
#include <fstream>
#include <algorithm>
#include <vector>
#include "../INPUT/Parameter.h"
#include "../MOLECULES/Molecules.h"
#include "util.h"


class Neighbors
{
    friend class Molecules;

private:
	std::vector<int> m_neighborList {};                // Neighbor list.
    std::vector<int> m_neighborIndex {};
    const double m_rSkin {};
    const double m_squareRSkin {};                        			// Skin radius squared.
    const int m_nDims {3};
    int m_updateRate {-1};
    const int m_numCell {};
    const double m_cellLength {};
	int m_errors { 0 };                                             // Errors of the neighbor list.
	std::vector<double> m_interDisplacementVector {};  // Inter neighbor list update displacement matrix.
    const std::vector<double> m_maxSquareRcArray{};
    const std::vector<double> m_threshArray{};
    int m_numNeighMax{};
    using NeighIterator = std::vector<int>::const_iterator;


public:
    // Molecule Neighbors. The neighbor list does NOT include bonded neighbors.
    Neighbors (param::Parameter param, const Molecules& systemMolecules)
            : m_rSkin { param.get_double( "rSkin") }
            , m_squareRSkin {std::pow (m_rSkin, 2 ) }
            , m_maxSquareRcArray (initializeMaxRc( systemMolecules))
            , m_threshArray (initializeThreshArray(param, m_maxSquareRcArray))
            , m_numCell { static_cast<int>(systemMolecules.m_lengthCube / m_rSkin) }
            , m_cellLength { systemMolecules.m_lengthCube / static_cast<double>(m_numCell)}

    {
        const double density {systemMolecules.m_nParticles / std::pow(systemMolecules.m_lengthCube, 3) };
        m_numNeighMax = static_cast<int>(std::pow(m_rSkin, 3) * 4./3. * 3.15 * density * 1.2);
        m_interDisplacementVector.resize(systemMolecules.m_nParticles * systemMolecules.m_nDims, 0);
        m_neighborList.resize(systemMolecules.m_nParticles);
        createNeighborList(systemMolecules);

        /***
        for(auto const& x : m_neighborList)
        {
            std::cout << x;
            if (it < m_numNeighMax) {
                it++;
                std::cout << " ";
            }
            else if (it == m_numNeighMax) {
                it=1;
                std::cout << "\n";
            }

        }
        ***/

    }


    static std::vector<double> initializeMaxRc(const Molecules& systemMolecules)
    {
        int nParticleTypes { systemMolecules.m_systemPairPotentials.getParticleTypes() };
        std::vector<double> maxSquareRcArray(nParticleTypes, 0) ;

        for (int i=1; i<=nParticleTypes; i++)
        {
            for (int j=i; j<=nParticleTypes; j++)
            {
                double squareRcIJ { systemMolecules.m_systemPairPotentials.getSquareRcIJ(i, j)};

                if (squareRcIJ > maxSquareRcArray[i-1])
                {
                    maxSquareRcArray[i-1] = squareRcIJ;
                }
                if (squareRcIJ > maxSquareRcArray[j-1])
                {
                    maxSquareRcArray[j-1] = squareRcIJ;
                }

            }
        }
        return maxSquareRcArray;
    };

    static std::vector<double> initializeThreshArray(param::Parameter param, const std::vector<double>& maxSquareRcArray)
    {
        double rSkin {param.get_double("rSkin")};
        std::vector<double> threshArray ;

        for (auto const &elt: maxSquareRcArray)
        {
            double maxRc { pow(elt, 1./2.) };
            const double rootThresh { (rSkin - maxRc)  / 2.};
            threshArray.push_back(rootThresh * rootThresh);
        }
        return threshArray;
    };

	void checkInterDisplacement(const Molecules& systemMolecules);


    [[nodiscard]] int cellTest(int indexCell) const;

   // void WOWcreateNeighborList(const Molecules& systemMolecules);

    [[nodiscard]] std::vector<std::vector<std::vector<std::vector<int>>>>
    createCellList(const Molecules &systemMolecules) const;

    void createSamePairCellNeighbor(const Molecules& systemMolecules,
                                    const std::vector<int>& oldNeighborList,
                                    const std::vector<int>& oldNeighborIndex,
                                    const std::vector<int>& cellMolecules, const bool& checkNeigh);

    void createDiffPairCellNeighbor(const Molecules& systemMolecules,
                                    const std::vector<int>& oldNeighborList,
                                    const std::vector<int>& oldNeighborIndex,
                                    const std::vector<int>& cellMolecules,
                                    const std::vector<int>& testNeighbors, const bool& checkNeigh);

    void checkNeighborError(const Molecules& systemMolecules, const std::vector<int> &oldNeighborList,
                            const std::vector<int>& oldNeighborIndex,
                            const double &squareDistance, const int &indexI, const int &indexJ);




    template<typename InputIt>
    void updateIJNeighbor(const Molecules& systemMolecules, const std::vector<int>& oldNeighborList,
                          const std::vector<int>& oldNeighborIndex, const InputIt& posItBegin,
                          const int& indexI, const int& indexJ, const bool& checkNeigh)
    {
        auto posItBeginJ {systemMolecules.getPosItBeginI(indexJ)};
        const double squareDistance { systemMolecules.squareDistancePair(posItBegin,
                                                                         posItBeginJ)};

        if (squareDistance < m_squareRSkin)
        {
            // j and i are neighbors
            const int lastIndexI {m_neighborIndex[indexI]};

            m_neighborList[lastIndexI] = indexJ; // j is on row i
            m_neighborIndex[indexI]++;

            const int lastIndexJ {m_neighborIndex[indexJ]};
            m_neighborList[lastIndexJ] = indexI; // j is on row i
            m_neighborIndex[indexJ]++;
            if (checkNeigh)
            {
                 checkNeighborError(systemMolecules, oldNeighborList, oldNeighborIndex,
                                    squareDistance, indexI, indexJ);
            }
        }
    }
    [[nodiscard]] int getNeighborIndexBegin(const int& particleIndex) const
    {
        //return m_neighborIndex[particleIndex];
        // std::cout  << particleIndex * m_numNeighMax <<'\n';
        return particleIndex * m_numNeighMax;
    }

    [[nodiscard]] int getNeighborIndexEnd(const int& particleIndex) const
    {
        return m_neighborIndex[particleIndex];
    }

    void sortNeighborList(const Molecules& systemMolecules);

    void createNeighborList(const Molecules& systemMolecules);

    void createCellNeighbors(const Molecules &systemMolecules,
                             const std::vector<int> &oldNeighborList,
                             const std::vector<int>& oldNeighborIndex,
                             const std::vector<std::vector<std::vector<std::vector<int>>>> &cellList, int xCell,
                             int yCell, int zCell, const bool &checkNeigh);

    [[nodiscard]] int getUpdateRate() const;

    [[nodiscard]] int getErrors() const;


    template<typename InputIt>
    void updateInterDisplacement(const int& indexTranslation, const InputIt& vectorItTranslation)
    {
        const int realIndex {indexTranslation * m_nDims};
        auto dispItBegin( m_interDisplacementVector.begin() + realIndex);

        std::transform(dispItBegin, dispItBegin + m_nDims, vectorItTranslation,
                       dispItBegin, std::plus<>());
    }

    [[nodiscard]] int getLenIndexBegin(const int &indexTranslation) const;

    [[nodiscard]] NeighIterator getNeighItBeginI(const int &indexParticle) const
    {
        const int& neighIndex {getNeighborIndexBegin(indexParticle)};
        return m_neighborList.begin() + neighIndex;
    }

    [[nodiscard]] NeighIterator getNeighItEndI(const int &indexParticle) const
    {
        const int& neighIndex {getNeighborIndexBegin(indexParticle + 1)};
        return m_neighborList.begin() + neighIndex;
    }
};


#endif /* NEIGHBORS_H_ */
