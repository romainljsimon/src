/*
 * MonteCarlo.cpp
 *
 *  Created on: 27 oct. 2022
 *      Author: Romain Simon
 */

#include <vector>
#include <cmath>
#include <algorithm>
#include "Neighbors.h"


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

void Neighbors::WOWcreateNeighborList(const Particles& systemParticles, const BondPotentials& systemBondPotentials)
{
    // Be careful with the swaps and poly dispersity. Solution for now is to take rSkin big enough.
    ++m_updateRate;
    std::vector<std::vector<int>> oldNeighborList = m_neighborList;
    m_neighborList.clear();
    m_neighborList.resize(m_nParticles);

    for (int i = 0; i < m_nParticles - 1; i++)
    {
        std::vector<double> positionParticle = systemParticles.getPositionI(i);

        for (int j = i + 1; j < m_nParticles; j++)
        {

            std::vector<int> bondsI {systemBondPotentials.getBondsI(i)};
            if (std::binary_search(bondsI.begin(), bondsI.end(), j))
            {
                continue;
            }

            double squareDistance {systemParticles.squareDistancePair(positionParticle, j)};

            if (squareDistance < m_squareRSkin)
            {
                // j and i are neighbors

                m_neighborList[i].push_back(j); // j is on row i
                m_neighborList[j].push_back(i); // i is on row j

                if (static_cast<int>(oldNeighborList[i].size()) > 0) // If this size is 0, then it is the first time the neighbor list is calculated.
                {

                    // Compare the new neighbor list with the new neighbor list

                    if (!std::binary_search(oldNeighborList[i].begin(), oldNeighborList[i].end(), j))
                    {
                        int typeI {systemParticles.getParticleTypeI(i) -1};
                        if (squareDistance <  m_maxSquareRcArray[typeI])
                        {
                            ++m_errors;
                        }
                    }
                }

            }
        }

    }
}

void Neighbors::createNeighborList(const Particles& systemParticles)
{
    // Be careful with the swaps and poly dispersity. Solution for now is to take rSkin big enough.
    const std::vector<std::vector<int>> oldNeighborList = m_neighborList;
    ++m_updateRate;
    m_neighborList.clear();
    m_neighborList.resize(m_nParticles);
    const std::vector<std::vector<std::vector<std::vector<int>>>> cellList { createCellList ( systemParticles) };
    const bool checkNeigh { m_updateRate > 0};

    for (int xCell = 0; xCell < m_numCell; xCell++)
    {
        for (int yCell = 0; yCell < m_numCell; yCell++)
        {
            for (int zCell = 0; zCell < m_numCell; zCell++)
            {
                createCellNeighbors(systemParticles, oldNeighborList, cellList, xCell, yCell, zCell, checkNeigh);
            }
        }
    }
    sortNeighborList();
}

void Neighbors::createNeighborList(const Particles& systemParticles, const BondPotentials& systemBondPotentials)
{
    // Be careful with the swaps and poly dispersity. Solution for now is to take rSkin big enough.
    const std::vector<std::vector<int>> oldNeighborList = m_neighborList;
    ++m_updateRate;
    m_neighborList.clear();
    m_neighborList.resize(m_nParticles);
    const std::vector<std::vector<std::vector<std::vector<int>>>> cellList { createCellList ( systemParticles) };
    const bool checkNeigh { m_updateRate > 0};

    for (int xCell = 0; xCell < m_numCell; xCell++)
    {
        for (int yCell = 0; yCell < m_numCell; yCell++)
        {
            for (int zCell = 0; zCell < m_numCell; zCell++)
            {
                createCellNeighbors(systemParticles, systemBondPotentials, oldNeighborList, cellList,
                                    xCell, yCell, zCell, checkNeigh);

            }
        }
    }
    sortNeighborList();
}

void Neighbors::sortNeighborList()
{
    for (int i = 0; i < m_nParticles; i++)
    {
        std::vector<int> neighborIList{m_neighborList[i]};
        std::sort(neighborIList.begin(), neighborIList.end());
        m_neighborList[i] = neighborIList;
    }
}

std::vector<std::vector<std::vector<std::vector<int>>>> Neighbors::createCellList(const Particles& systemParticles) const
{
    std::vector<std::vector<std::vector<std::vector<int>>>> cellList{};
    cellList.resize(m_numCell, std::vector<std::vector<std::vector<int>>>(m_numCell,
                                                                          std::vector<std::vector<int>>(m_numCell)));

    for (int i = 0; i < m_nParticles; i++)
    {
        std::vector<double> positionArrayI { systemParticles.getPositionI(i) };
        const int xCell{static_cast<int>(floor(positionArrayI[0] / m_cellLength))};
        const int yCell{static_cast<int>(floor(positionArrayI[1] / m_cellLength))};
        const int zCell{static_cast<int>(floor(positionArrayI[2] / m_cellLength))};
        cellList[xCell][yCell][zCell].push_back(i);
    }

    return cellList;
}

void Neighbors::createCellNeighbors(const Particles& systemParticles,
                                    const std::vector<std::vector<int>>& oldNeighborList,
                                    const std::vector<std::vector<std::vector<std::vector<int>>>>& cellList,
                                    int xCell, int yCell, int zCell, const bool& checkNeigh)
{
    const std::vector<int> xyzCellParticles{cellList[xCell][yCell][zCell]};
    createSamePairCellNeighbor(systemParticles, oldNeighborList, xyzCellParticles, checkNeigh);
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
                    //std::cout << xCell << "  " << testXCell << "\n";
                    // std::cout << yCell << "  " << testYCell << "\n";
                    // std::cout << zCell << "  " << testZCell << "\n";
                    createDiffPairCellNeighbor(systemParticles,
                                               oldNeighborList, xyzCellParticles,
                                               cellList[testXCell][testYCell][testZCell], checkNeigh);
                }
            }
        }
    }
}

void Neighbors::createCellNeighbors(const Particles& systemParticles, const BondPotentials& systemBondPotentials,
                                    const std::vector<std::vector<int>>& oldNeighborList,
                                    const std::vector<std::vector<std::vector<std::vector<int>>>>& cellList,
                                    int xCell, int yCell, int zCell, const bool& checkNeigh)
{
    const std::vector<int> xyzCellParticles{cellList[xCell][yCell][zCell]};
    createSamePairCellNeighbor(systemParticles, systemBondPotentials, oldNeighborList,
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
                    //std::cout << xCell << "  " << testXCell << "\n";
                    // std::cout << yCell << "  " << testYCell << "\n";
                    // std::cout << zCell << "  " << testZCell << "\n";
                    createDiffPairCellNeighbor(systemParticles, systemBondPotentials,
                                               oldNeighborList, xyzCellParticles,
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


void Neighbors::createSamePairCellNeighbor(const Particles& systemParticles,
                                           const std::vector<std::vector<int>>& oldNeighborList,
                                           const std::vector<int>& cellParticles, const bool& checkNeigh)
{
    const int cellParticlesSize{static_cast<int>(cellParticles.size())};

    for (int i = 0; i < cellParticlesSize - 1; i++)
    {
        const int realIIndex { cellParticles[i] };
        const std::vector<double> positionParticle { systemParticles.getPositionI(realIIndex) };

        for (int j = i+1; j < cellParticlesSize; j++)
        {
            const int realJIndex {cellParticles[j]};
            updateIJNeighbor(systemParticles, oldNeighborList, positionParticle,
                             realIIndex, realJIndex, checkNeigh);

        }
    }
}

void Neighbors::createSamePairCellNeighbor(const Particles& systemParticles, const BondPotentials& systemBondPotentials,
                                           const std::vector<std::vector<int>>& oldNeighborList,
                                           const std::vector<int>& cellParticles, const bool& checkNeigh)
{
    const int cellParticlesSize{static_cast<int>(cellParticles.size())};

    for (int i = 0; i < cellParticlesSize - 1; i++)
    {
        const int realIIndex { cellParticles[i] };

        const std::vector<double> positionParticle { systemParticles.getPositionI(realIIndex) };

        const std::vector<int> bondsParticleI{ systemBondPotentials.getBondsI(realIIndex) };

        for (int j = i+1; j < cellParticlesSize; j++)
        {
            const int realJIndex {cellParticles[j]};


            if (!std::binary_search(bondsParticleI.begin(), bondsParticleI.end(), realJIndex))
            {
                updateIJNeighbor(systemParticles, oldNeighborList, positionParticle,
                                 realIIndex, realJIndex, checkNeigh);
            }
        }
    }
}

void Neighbors::createDiffPairCellNeighbor(const Particles& systemParticles,
                                           const std::vector<std::vector<int>>& oldNeighborList,
                                           const std::vector<int>& cellParticles,
                                           const std::vector<int>& testNeighbors, const bool& checkNeigh)
{
    for (auto const &i: cellParticles)
    {
        const std::vector<double> positionParticle { systemParticles.getPositionI(i) };

        for (auto const &j: testNeighbors)
        {
            updateIJNeighbor(systemParticles, oldNeighborList, positionParticle, i, j, checkNeigh);
        }
    }
}

void Neighbors::createDiffPairCellNeighbor(const Particles& systemParticles, const BondPotentials& systemBondPotentials,
                                           const std::vector<std::vector<int>>& oldNeighborList,
                                           const std::vector<int>& cellParticles,
                                           const std::vector<int>& testNeighbors, const bool& checkNeigh)
{
    for (auto const &i: cellParticles)
    {
        const std::vector<double> positionParticle { systemParticles.getPositionI(i) };
        const std::vector<int> bondsParticleI{ systemBondPotentials.getBondsI(i) };

        for (auto const &j: testNeighbors)
        {
            if (!std::binary_search(bondsParticleI.begin(), bondsParticleI.end(), j))
            {
                updateIJNeighbor(systemParticles, oldNeighborList, positionParticle, i, j, checkNeigh);
            }
        }
    }
}

void Neighbors::updateIJNeighbor(const Particles& systemParticles, const std::vector<std::vector<int>>& oldNeighborList,
                                 const std::vector<double>& positionParticle, const int& indexI, const int& indexJ,
                                 const bool& checkNeigh)
{
    const double squareDistance { systemParticles.squareDistancePair(positionParticle, indexJ) };

    if (squareDistance < m_squareRSkin)
    {
        // j and i are neighbors
        m_neighborList[indexI].push_back(indexJ); // j is on row i
        m_neighborList[indexJ].push_back(indexI); // i is on row j
        if (checkNeigh)
        {
            checkNeighborError(systemParticles, oldNeighborList, squareDistance, indexI, indexJ);
        }
    }
}



void Neighbors::checkNeighborError(const Particles& systemParticles, const std::vector<std::vector<int>>& oldNeighborList,
                                   const double& squareDistance, const int& indexI, const int& indexJ)
{
    //if (static_cast<int>(oldNeighborList[indexI].size()) > 0) // If this size is 0, then it is the first time the neighbor list is calculated.
    //{
        // Compare the new neighbor list with the new neighbor list
    if (!std::binary_search(oldNeighborList[indexI].begin(), oldNeighborList[indexI].end(), indexJ))
    {

        int particleTypeIndex { systemParticles.getParticleTypeI(indexI) - 1};

        // double max_rc{ (m_maxDiam + m_typeArray[indexI]) / 2 * m_rC};
        double squareMaxRc { m_maxSquareRcArray[particleTypeIndex] };
        if (squareDistance < squareMaxRc)
        {
            ++m_errors;
        }
    }
    //}
}



void Neighbors::updateInterDisplacement(const int& indexTranslation, const std::vector<double>& vectorTranslation)
{
    m_interDisplacementMatrix[indexTranslation] = vectorSum(m_interDisplacementMatrix[indexTranslation],
                                                            vectorTranslation);
}

/*******************************************************************************
 * This function checks if the maximum displacement of one particle is superior
 * than m_squareRDiff. If this is true, then the neighbor list is updated and
 * m_interDisplacementMatrix is reinitialized to zero. The criterion shouldn't
 * allow for any errors?
 ******************************************************************************/
void Neighbors::checkInterDisplacement(const Particles& systemParticles)
{
    const  std::vector<double> squareDispVector = getSquareNormRowMatrix(m_interDisplacementMatrix);
    int nParticles { systemParticles.getNParticles()};
    for (int i=0; i < nParticles; i++)
    {
        int particleTypeIndex { systemParticles.getParticleTypeI(i) - 1};
        //const double maxRc {( m_maxDiam + m_typeArray[i]) / 2 * m_rC};
        //const double thresh = { pow ( (m_rSkin - maxRc) / 2, 2)};
        if ( squareDispVector[i] > m_threshArray[particleTypeIndex])
        {
            createNeighborList(systemParticles);
            std::fill(m_interDisplacementMatrix.begin(), m_interDisplacementMatrix.end(),
                      std::vector<double>(3, 0));
            break;
        }
    }
}

void Neighbors::checkInterDisplacement(const Particles& systemParticles, const BondPotentials& systemBondPotentials)
{
    const std::vector<double> squareDispVector { getSquareNormRowMatrix(m_interDisplacementMatrix) };
    int nParticles { systemParticles.getNParticles()};
    for (int i=0; i < nParticles; i++)
    {
        int particleTypeIndex { systemParticles.getParticleTypeI(i) - 1};
        //const double maxRc {( m_maxDiam + m_typeArray[i]) / 2 * m_rC};
        //const double thresh = { pow ( (m_rSkin - maxRc) / 2, 2)};
        const double& thresh {m_threshArray[particleTypeIndex]};
        if ( squareDispVector[i] > thresh)
        {
            WOWcreateNeighborList(systemParticles, systemBondPotentials);
            std::fill(m_interDisplacementMatrix.begin(), m_interDisplacementMatrix.end(),
                      std::vector<double>(3, 0));
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

std::vector<int> Neighbors::getNeighborIList(int particleIndex) const
{
    /***
    std::vector<int> neighborIList { };

    if ( m_neighMethod == "verlet" )
    {neighborIList = m_neighborList[particleIndex]; // particle's neighbor list
    }

    else
    {
        neighborIList.resize ( m_nParticles );
        std::iota (std::begin(neighborIList), std::end(neighborIList), 0);
    }
    ***/
    return m_neighborList[particleIndex];
}

