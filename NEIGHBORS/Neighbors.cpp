/*
 * MonteCarlo.cpp
 *
 *  Created on: 27 oct. 2022
 *      Author: Romain Simon
 */

#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include "Neighbors.h"
#include "../MOLECULES/Molecules.h"


/*******************************************************************************
* VERLET LIST AND CELL LIST CREATION AND UPDATE METHODS
 ******************************************************************************/


/*******************************************************************************
 * Creation or update of the Verlet neighbor list. The neighbor list is
 * implemented as an 2d array such as row i of the array is particle i's
 * neighbor list. If j is on row i then i is on row j.
 *
 * Example of a neighbor list with 4 particles:
 * row
 *  1: 2 3
 *  2: 1 3 4
 *  3: 1 2
 *  4: 2
 *
 * Each time the neighbor list is updated, it is compared with the previous
 * neighbor list to detect potential errors. Normally, the code has been made
 * such as there are no errors.
 ******************************************************************************/
/***
void Neighbors::WOWcreateNeighborList(const Molecules& systemMolecules)
{
    // Be careful with the swaps and poly dispersity. Solution for now is to take rSkin big enough.
    ++m_updateRate;
    const std::vector<std::vector<int>> oldNeighborList { m_neighborList};
    m_neighborList.clear();
    m_neighborList.resize(systemMolecules.m_nParticles);

    for (int i = 0; i < systemMolecules.m_nParticles - 1; i++)
    {
        const std::vector<double>& positionParticle { systemMolecules.getPositionI(i) };
        const std::vector<int>& bondsI {systemMolecules.m_bondsArray[i]};
        //const int bondsSize {static_cast<int>(bondsI.size())};
        //int count {0};

        for (int j = i + 1; j < systemMolecules.m_nParticles; j++)
        {

            const bool notIn {!std::binary_search(bondsI.begin(), bondsI.end(), j)};
            if (notIn)
            {


                const double squareDistance{systemMolecules.squareDistancePair(systemMolecules.m_positionArray.begin()
                                                                               + 3 * i,
                                                                               systemMolecules.m_positionArray.begin()
                                                                               + 3 * j)};

                if (squareDistance < m_squareRSkin) {
                    // j and i are neighbors

                    m_neighborList[i].push_back(j); // j is on row i
                    m_neighborList[j].push_back(i); // i is on row j

                    if (static_cast<int>(oldNeighborList[i].size()) > 0) // If this size is 0, then it is the first time the neighbor list is calculated.
                    {

                        // Compare the new neighbor list with the old neighbor list

                        if (!std::binary_search(oldNeighborList[i].begin(), oldNeighborList[i].end(), j)) {
                            const int typeI{systemMolecules.getParticleTypeI(i) - 1};
                            if (squareDistance < m_maxSquareRcArray[typeI])
                            {
                                ++m_errors;
                            }
                        }
                    }

                }
            }
        }

    }
}
***/
void Neighbors::createNeighborList(const Molecules& systemMolecules)
{
    // Be careful with the swaps and poly dispersity. Solution for now is to take rSkin big enough.
    const std::vector<int> oldNeighborList = m_neighborList;
    const std::vector<int> oldNeighborIndex = m_neighborIndex;
    ++m_updateRate;

    m_neighborList.clear();
    m_neighborList.resize(systemMolecules.m_nParticles * m_numNeighMax, -1);
    m_neighborIndex.clear();
    m_neighborIndex.resize(systemMolecules.m_nParticles);
    std::iota(m_neighborIndex.begin(), m_neighborIndex.end(), 0);
    std::for_each(m_neighborIndex.begin(), m_neighborIndex.end(), [this](int &n)
    { n=n*m_numNeighMax; });
    // std::cout << m_neighborIndex[1] <<"\n";
    //for (int i = 0; i<systemMolecules.m_nParticles; i++)
    //{
    //    m_neighborList[i].clear();
    //}
    const std::vector<std::vector<std::vector<std::vector<int>>>> cellList { createCellList ( systemMolecules) };
    const bool checkNeigh { m_updateRate > 0};

    for (int xCell = 0; xCell < m_numCell; xCell++)
    {
        for (int yCell = 0; yCell < m_numCell; yCell++)
        {
            for (int zCell = 0; zCell < m_numCell; zCell++)
            {
                createCellNeighbors(systemMolecules, oldNeighborList, oldNeighborIndex,
                                    cellList, xCell, yCell, zCell, checkNeigh);
            }
        }
    }
    // sortNeighborList(systemMolecules);
}

void Neighbors::sortNeighborList(const Molecules& systemMolecules)
{
    for (int i = 0; i < systemMolecules.m_nParticles; i++)
    {
        std::sort(m_neighborList.begin() + i * m_numNeighMax, m_neighborList.begin() + m_neighborIndex[i]);
    }
}

std::vector<std::vector<std::vector<std::vector<int>>>> Neighbors::createCellList(const Molecules& systemMolecules) const
{
    std::vector<std::vector<std::vector<std::vector<int>>>> cellList{};
    cellList.resize(m_numCell, std::vector<std::vector<std::vector<int>>>(m_numCell,
                                                                          std::vector<std::vector<int>>(m_numCell)));

    for (int i = 0; i < systemMolecules.m_nParticles; i++)
    {
        auto posItBeginI { systemMolecules.getPosItBeginI(i)};
        const int xCell{static_cast<int>(floor(*posItBeginI / m_cellLength))};
        posItBeginI++;
        const int yCell{static_cast<int>(floor(*posItBeginI / m_cellLength))};
        posItBeginI++;
        const int zCell{static_cast<int>(floor(*posItBeginI / m_cellLength))};
        cellList[xCell][yCell][zCell].push_back(i);
    }

    return cellList;
}


void Neighbors::createCellNeighbors(const Molecules& systemMolecules,
                                    const std::vector<int>& oldNeighborList,
                                    const std::vector<int>& oldNeighborIndex,
                                    const std::vector<std::vector<std::vector<std::vector<int>>>>& cellList,
                                    int xCell, int yCell, int zCell, const bool& checkNeigh)
{
    const std::vector<int> xyzCellParticles{ cellList[xCell][yCell][zCell] };
    createSamePairCellNeighbor(systemMolecules, oldNeighborList, oldNeighborIndex,
                               xyzCellParticles, checkNeigh);

    const int xyzInt {xCell * 100 + yCell * 10 + zCell};
    for (int xCellDiff = -1; xCellDiff < 2; xCellDiff++)
    {
        const int testXCell {cellTest(xCell + xCellDiff)};

        for (int yCellDiff = -1; yCellDiff < 2; yCellDiff++)
        {
            const int testYCell {cellTest(yCell + yCellDiff)};

            for (int zCellDiff = -1; zCellDiff < 2; zCellDiff++)
            {
                const int testZCell {cellTest(zCell + zCellDiff)};

                const int testInt {testXCell * 100 + testYCell * 10 + testZCell};

                if (testInt  > xyzInt)
                {

                    createDiffPairCellNeighbor(systemMolecules,
                                               oldNeighborList, oldNeighborIndex, xyzCellParticles,
                                               cellList[testXCell][testYCell][testZCell], checkNeigh);
                }
            }
        }
    }
}


int Neighbors::cellTest(int indexCell) const
{
    if (indexCell < 0)
    {
        return m_numCell - 1;
    }
    else if (indexCell == m_numCell)
    {
        return 0;
    }
    else
    {
        return indexCell;
    }
}

void Neighbors::createSamePairCellNeighbor(const Molecules& systemMolecules,
                                           const std::vector<int>& oldNeighborList,
                                           const std::vector<int>& oldNeighborIndex,
                                           const std::vector<int>& cellParticles, const bool& checkNeigh)
{
    const int cellParticlesSize{static_cast<int>(cellParticles.size())};

    for (int i = 0; i < cellParticlesSize - 1; i++)
    {
        const int realIIndex{cellParticles[i]};

        const auto &posItBegin{systemMolecules.getPosItBeginI(realIIndex)};

        //const std::vector<int> &bondsParticleI{systemMolecules.m_bondsArray[realIIndex]};
        const auto& bondsItBegin {systemMolecules.getBondsItBeginI(realIIndex)};
        const auto& bondsItEnd {systemMolecules.getBondsItEndI(realIIndex)};
        for (int j = i + 1; j < cellParticlesSize; j++) {
            const int realJIndex{cellParticles[j]};


            if (!std::binary_search(bondsItBegin, bondsItEnd, realJIndex))
            {
                updateIJNeighbor(systemMolecules, oldNeighborList, oldNeighborIndex, posItBegin,
                                 realIIndex, realJIndex, checkNeigh);
            }
        }
    }
}

void Neighbors::createDiffPairCellNeighbor(const Molecules& systemMolecules,
                                           const std::vector<int>& oldNeighborList,
                                           const std::vector<int>& oldNeighborIndex,
                                           const std::vector<int>& cellParticles,
                                           const std::vector<int>& testNeighbors, const bool& checkNeigh)
{
    for (auto const &i: cellParticles)
    {
        const auto& posItBegin { systemMolecules.getPosItBeginI(i)};
        const auto& bondsItBegin {systemMolecules.getBondsItBeginI(i)};
        const auto& bondsItEnd {systemMolecules.getBondsItEndI(i)};

        for (auto const &j: testNeighbors)
        {
            if (!std::binary_search(bondsItBegin, bondsItEnd, j))
            {
                updateIJNeighbor(systemMolecules, oldNeighborList, oldNeighborIndex, posItBegin, i, j, checkNeigh);
            }
        }
    }
}




void Neighbors::checkNeighborError(const Molecules& systemMolecules, const std::vector<int>& oldNeighborList,
                                   const std::vector<int>& oldNeighborIndex,
                                   const double& squareDistance, const int& indexI, const int& indexJ)
{
    //if (static_cast<int>(oldNeighborList[indexI].size()) > 0) // If this size is 0, then it is the first time the neighbor list is calculated.
    //{
        // Compare the new neighbor list with the new neighbor list
    const auto& oldItBegin { oldNeighborList.begin() + indexI * m_numNeighMax };
    const auto& oldItEnd { oldNeighborList.begin() + oldNeighborIndex[indexI]};

    //if (!std::binary_search(oldItBegin, oldItEnd, indexJ))
    if ((std::find(oldItBegin, oldItEnd, indexJ) == oldItEnd))
    {
        int particleTypeIndex { systemMolecules.getParticleTypeI(indexI) - 1};

        double squareMaxRc { m_maxSquareRcArray[particleTypeIndex] };
        if (squareDistance < squareMaxRc)
        {
            ++m_errors;
        }
    }
    //}
}



/*******************************************************************************
 * This function checks if the maximum displacement of one particle is superior
 * than m_squareRDiff. If this is true, then the neighbor list is updated and
 * m_interDisplacementMatrix is reinitialized to zero. The criterion shouldn't
 * allow for any errors?
 ******************************************************************************/
void Neighbors::checkInterDisplacement(const Molecules& systemMolecules)
{
    const  std::vector<double> squareDispVector = getSquareNormRowMatrix(m_interDisplacementVector.begin(),
                                                                         systemMolecules.m_nParticles,
                                                                         m_nDims);


    for (int i=0; i < systemMolecules.m_nParticles; i++)
    {
        int particleTypeIndex { systemMolecules.getParticleTypeI(i) - 1};

        if ( squareDispVector[i] > m_threshArray[particleTypeIndex])
        {
            createNeighborList(systemMolecules);
            std::fill(m_interDisplacementVector.begin(), m_interDisplacementVector.end(),0);
            break;
        }
    }
}


/*******************************************************************************
* EXTRACT NEIGHBOR INFORMATION METHODS
 ******************************************************************************/

int Neighbors::getUpdateRate() const
{
    return m_updateRate;
}

int Neighbors::getErrors() const
{
    return m_errors;
}

int Neighbors::getLenIndexBegin(const int& indexParticle) const
{
    const int lenNeigh {getNeighborIndexEnd(indexParticle) - getNeighborIndexBegin(indexParticle)};
    return lenNeigh;
}


