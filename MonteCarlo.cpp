/*
 * MonteCarlo.cpp
 *
 *  Created on: 27 oct. 2022
 *      Author: Romain Simon
 */

#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <numeric>
#include <algorithm>
#include "MonteCarlo.h"
#include "random.h"
#include "energy.h"
#include "readSaveFile.h"
#include "util.h"


/*******************************************************************************
 * This function is the core of the Monte Carlo program. It iterates over
 * m_timeSteps steps that the user defines. At each step m_nParticles Monte
 * Carlo moves are made. For now the Monte Carlo move is a translation of a
 * randomly chosen particle.
 * The energy, particles positions, particles displacements and pressure
 * (optional) are written in files.
 ******************************************************************************/
void MonteCarlo::mcTotal()
{
	std::string preName ("./outXYZ/position");
	std::string extname {".xyz"};
	std::string preNameDisp ("./disp/displacement");
	std::string extnameDisp {".txt"};
    std::vector<double> radiusArray (divideVectorByScalar(m_diameterArray, 2));

	saveInXYZ(m_positionArray, radiusArray, m_moleculeType, m_lengthCube,preName + std::to_string(0) + extname );
	saveDoubleTXT(m_energy / m_nParticles, m_folderPath + "/outE.txt");
	saveDisplacement(m_totalDisplacementMatrix, preNameDisp + std::to_string(0) + extnameDisp);

	const std::vector<int> saveTimeStepArray ( createSaveTime(m_timeSteps, m_saveUpdate, 1.1));

	int save_index { 0 };


	for (int i = 0; i < m_timeSteps; i++) //Iteration over m_timeSteps
	{

		for (int j = 0; j < m_nParticles; j++) // N Monte Carlo moves are tried in one time step.
		{
			mcMove();
		}

		m_interDisplacementMatrix = matrixSum( m_interDisplacementMatrix, m_stepDisplacementMatrix );
		m_totalDisplacementMatrix = matrixSum( m_totalDisplacementMatrix, m_stepDisplacementMatrix );
		std::fill( m_stepDisplacementMatrix.begin(), m_stepDisplacementMatrix.end(), std::vector<double>(3, 0));

		if (m_neighMethod == "verlet")
		{
			checkStepDisplacement();
		}

		// Next, the results of the simulations are saved.

		if (saveTimeStepArray[save_index] == i)
		{
			radiusArray = divideVectorByScalar(m_diameterArray, 2);
            saveInXYZ (m_positionArray, radiusArray, m_moleculeType, m_lengthCube, preName + std::to_string(i + 1) + extname );
			saveDisplacement (m_totalDisplacementMatrix, preNameDisp + std::to_string(i + 1) + extnameDisp );
			++save_index;
		}

		if (i % m_saveRate == 0)
		{
			saveDoubleTXT(m_energy / m_nParticles, m_folderPath + "/outE.txt"); //Energy is saved at each time step.
		}
		if (m_calculatePressure) // Saves pressure if m_calculatePressure is True
		{
			saveDoubleTXT( m_pressure, m_folderPath + "/outP.txt");
		}

	}

    radiusArray = divideVectorByScalar(m_diameterArray, 2);
	m_acceptanceRate /= m_timeSteps;
    saveInXYZ(m_positionArray, radiusArray, m_moleculeType, m_lengthCube, preName + std::to_string(m_timeSteps) + extname );
    saveDisplacement(m_totalDisplacementMatrix, preNameDisp + std::to_string(m_timeSteps) + extnameDisp);
    saveDoubleTXT( m_errors, m_folderPath + "/errors.txt");
    std::cout << "Translation MC move acceptance rate: " << m_acceptanceRate - m_pSwap * m_acceptanceRate << "\n";

    if ( m_swap )
    {
        m_acceptanceRateSwap /= m_timeSteps;
        m_acceptanceRateSwap /= m_pSwap;
        std::cout << "Swap MC move acceptance rate: " << m_acceptanceRateSwap << "\n";
        std::cout << "Total MC move acceptance rate: " << m_acceptanceRate << "\n";
    }

	std::cout << "Neighbor list update rate: " << m_updateRate / m_timeSteps << "\n";
	std::cout << "Number of neighbor list errors: " << m_errors << "\n";

}


std::vector<std::vector<std::vector<std::vector<int>>>> MonteCarlo::createCellList()
{
    std::vector<std::vector<std::vector<std::vector<int>>>> cellList{};
    cellList.resize(m_numCell, std::vector<std::vector<std::vector<int>>>(m_numCell, std::vector<std::vector<int>>(m_numCell)));

    for (int i = 0; i < m_nParticles; i++)
    {
        const int xCell{static_cast<int>(floor(m_positionArray[i][0] / m_cellLength))};
        const int yCell{static_cast<int>(floor(m_positionArray[i][1] / m_cellLength))};
        const int zCell{static_cast<int>(floor(m_positionArray[i][2] / m_cellLength))};
        cellList[xCell][yCell][zCell].push_back(i);
    }

    return cellList;
}

int MonteCarlo::cellTest(int indexCell) const
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
void MonteCarlo::checkNeighborError(const std::vector<std::vector<int>>& oldNeighborList, const double& squareDistance, const int& indexI, const int& indexJ)
{
    if (static_cast<int>(oldNeighborList[indexI].size()) > 0) // If this size is 0, then it is the first time the neighbor list is calculated.
    {
        // Compare the new neighbor list with the new neighbor list
        if (!std::binary_search(oldNeighborList[indexI].begin(), oldNeighborList[indexI].end(), indexJ))
        {
            double max_rc{(m_maxDiam + m_diameterArray[indexI]) / 2 * m_rC};

            if (squareDistance < pow(max_rc, 2))
            {
                ++m_errors;
            }
        }
    }
}
void MonteCarlo::updateIJNeighbor(const std::vector<std::vector<int>>& oldNeighborList, const std::vector<double>& positionParticle, const int& indexI, const int& indexJ, const bool& checkNeigh)
{
    const double squareDistance{
            squareDistancePair(positionParticle, m_positionArray[indexJ], m_lengthCube, m_halfLengthCube)};

    if (squareDistance < m_squareRSkin)
    {
        // j and i are neighbors
        m_neighborList[indexI].push_back(indexJ); // j is on row i
        m_neighborList[indexJ].push_back(indexI); // i is on row j
        if (checkNeigh)
        {
            checkNeighborError(oldNeighborList, squareDistance, indexI, indexJ);
        }
    }
}

void MonteCarlo::createSamePairCellNeighbor(const std::vector<std::vector<int>>& oldNeighborList,
                                            const std::vector<int>& cellParticles, const bool& checkNeigh)
{
    const int cellParticlesSize{static_cast<int>(cellParticles.size())};

    for (int i = 0; i < cellParticlesSize - 1; i++)
    {
        const int realIIndex {cellParticles[i]};
        const std::vector<double> positionParticle{m_positionArray[realIIndex]};

        for (int j = i+1; j < cellParticlesSize; j++)
        {
            const int realJIndex {cellParticles[j]};
            updateIJNeighbor(oldNeighborList, positionParticle, realIIndex, realJIndex, checkNeigh);

        }
    }
}

void MonteCarlo::createDiffPairCellNeighbor(const std::vector<std::vector<int>>& oldNeighborList,
                                            const std::vector<int>& cellParticles,
                                            const std::vector<int>& testNeighbors, const bool& checkNeigh)
{
    for (auto const &i: cellParticles)
    {
        const std::vector<double> positionParticle{m_positionArray[i]};

        for (auto const &j: testNeighbors)
        {
            updateIJNeighbor(oldNeighborList, positionParticle, i, j, checkNeigh);
        }
    }
}

void MonteCarlo::createCellNeighbors(const std::vector<std::vector<int>>& oldNeighborList,
                                     const std::vector<std::vector<std::vector<std::vector<int>>>>& cellList,
                                     int xCell, int yCell, int zCell, const bool& checkNeigh)
{
    const std::vector<int> xyzCellParticles{cellList[xCell][yCell][zCell]};
    createSamePairCellNeighbor(oldNeighborList, xyzCellParticles, checkNeigh);
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
                    createDiffPairCellNeighbor(oldNeighborList, xyzCellParticles,
                                               cellList[testXCell][testYCell][testZCell], checkNeigh);
                }


            }
        }


    }
}
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
void MonteCarlo::createNeighborList()
{
    // Be careful with the swaps and poly dispersity. Solution for now is to take rSkin big enough.
	++m_updateRate;
	const std::vector<std::vector<int>> oldNeighborList = m_neighborList;
	m_neighborList.clear();
	m_neighborList.resize(m_nParticles);
    const std::vector<std::vector<std::vector<std::vector<int>>>> cellList { createCellList() };
    const bool checkNeigh { m_updateRate > 0};
	for (int xCell = 0; xCell < m_numCell; xCell++)
	{
        for (int yCell = 0; yCell < m_numCell; yCell++)
        {
            for (int zCell = 0; zCell < m_numCell; zCell++)
            {
                createCellNeighbors(oldNeighborList, cellList, xCell, yCell, zCell, checkNeigh);
            }
		}
	}
    for (int i = 0; i < m_nParticles; i++)
    {
        std::vector<int> neighborIList{m_neighborList[i]};
        std::sort(neighborIList.begin(), neighborIList.end());
        m_neighborList[i] = neighborIList;
    }
}

/*******************************************************************************
 * This function implements a Monte Carlo move: translation of a random particle,
 * calculation of the energy of the new system and then acceptation or not of
 * the move according to the Metropolis criterion. It is called N times in one
 * time step.
 ******************************************************************************/
void MonteCarlo::mcMove()
{
    if ( m_swap )
    {
        const double randomDouble { randomDoubleGenerator(0., 1.) } ;
        const bool swapped {randomDouble < m_pSwap};

        if ( swapped )
        {
            mcSwap();
        }
        else
        {
            mcTranslation ();
        }
    }
    else
    {
        mcTranslation ();
    }
}

/*******************************************************************************
 * This function implements a Monte Carlo move: translation of a the whole
 * polymer, calculation of the energy of the new system and then acceptation or not of
 * the move according to the Metropolis criterion. It is called ??? times in one
 * time step.
 ******************************************************************************/
/***
void MonteCarlo::mcAllMove()
{
	int indexTranslation { randomIntGenerator(0, (m_nParticles - 1) / 3) }; // randomly chosen particle
	std::vector<double> randomVector ( randomVectorDoubleGenerator(3, -m_rBox, m_rBox) );

	double oldEnergyPolymer {};
	double newEnergyPolymer {};
	std::vector< std::vector<double>> tentativePositionArray(m_positionArray);

	for (int j = 0; j < 2; j++)
	{

		std::vector<double> positionParticleTranslation = mcTranslation ( indexTranslation + j, randomVector );

		tentativePositionArray[indexTranslation + j] = positionParticleTranslation;

	}

	for (int j = 0; j < 2; j++)
	{
		std::vector<int> neighborIList { };

		if ( m_neighMethod == "verlet" )
		{
			neighborIList = m_neighborList[indexTranslation + j]; // particle's neighbor list
		}

		else
		{
			neighborIList.resize ( m_nParticles );
			std::iota (std::begin(neighborIList), std::end(neighborIList), 0);

		}

		double oldEnergyParticle {};
		double newEnergyParticle {};

 		std::vector<int> bondsI ( m_bondsMatrix[indexTranslation] );

		oldEnergyParticle = energyParticlePolymer (indexTranslation, m_positionArray[indexTranslation],
												   m_positionArray, neighborIList, m_diameterArray, bondsI,
												   m_squareRc, m_lengthCube, m_squareR0, m_feneK);
		newEnergyParticle = energyParticlePolymer (indexTranslation, tentativePositionArray[indexTranslation + j],
												   tentativePositionArray, neighborIList, m_diameterArray, bondsI,
												   m_squareRc, m_lengthCube, m_squareR0, m_feneK);

		oldEnergyPolymer += oldEnergyParticle;
		newEnergyParticle += newEnergyParticle;

	}


	// Metropolis criterion
	bool acceptMove { metropolis (newEnergyPolymer, oldEnergyPolymer) } ;

	// If the move is accepted, then the energy, the position array and the displacement array can be updated.
	// If the m_calculatePressure is set to True, then the pressure is calculated.
	if ( acceptMove )
	{



		double newEnergy {m_energy - oldEnergyPolymer + newEnergyPolymer};
		m_energy = newEnergy;
		m_stepDisplacementMatrix[indexTranslation] = vectorSum(m_stepDisplacementMatrix[indexTranslation], randomVector);
		m_positionArray = tentativePositionArray;
		// m_acceptanceRate += 1. / m_nParticles; // increment of the acceptance rate.

	}

}
***/
/*******************************************************************************
 * This function returns a tentative new particle position.
 *
 * @param indexTranslation Translated particle's index.
 *        randomVector Random vector of size 3. The components of the vector are
 *                     taken from a uniform distribution U(-m_rBox, m_rBox).
 * @return Tentative new particle position.
 ******************************************************************************/
void MonteCarlo::mcTranslation() {
    const int indexTranslation{randomIntGenerator(0, m_nParticles - 1)}; // randomly chosen particle
    const std::vector<double> randomVector(randomVectorDoubleGenerator(3, -m_rBox, m_rBox));
    std::vector<double> positionParticleTranslation = m_positionArray[indexTranslation];
    positionParticleTranslation = vectorSum(positionParticleTranslation, randomVector);
    positionParticleTranslation = periodicBC(positionParticleTranslation, m_lengthCube);

    const std::vector<int> neighborIList{findNeighborIList( indexTranslation )};

    double oldEnergyParticle;
    double newEnergyParticle;

    if (m_simulationMol == "polymer")
    {
        const std::vector<int> bondsI(m_bondsMatrix[indexTranslation]);

        oldEnergyParticle = energyParticlePolymer(indexTranslation, m_positionArray[indexTranslation],
                                                  m_positionArray, neighborIList, m_diameterArray, bondsI,
                                                  m_squareRc, m_lengthCube, m_halfLengthCube,m_squareR0, m_feneK, m_bondType);
        newEnergyParticle = energyParticlePolymer(indexTranslation, positionParticleTranslation,
                                                  m_positionArray, neighborIList, m_diameterArray, bondsI,
                                                  m_squareRc, m_lengthCube, m_halfLengthCube, m_squareR0, m_feneK, m_bondType);

    }
    else
    {
        oldEnergyParticle = energyParticle(indexTranslation, m_positionArray[indexTranslation], m_positionArray,
                                           neighborIList, m_diameterArray, m_squareRc, m_lengthCube, m_halfLengthCube);
        newEnergyParticle = energyParticle(indexTranslation, positionParticleTranslation, m_positionArray,
                                           neighborIList, m_diameterArray, m_squareRc, m_lengthCube, m_halfLengthCube);
    }
    const double diff_energy {newEnergyParticle - oldEnergyParticle};

    // Metropolis criterion
    const bool acceptMove{ metropolis(diff_energy) };

    // If the move is accepted, then the energy, the position array and the displacement array can be updated.
    // If the m_calculatePressure is set to True, then the pressure is calculated.
    if (acceptMove)
    {
        generalUpdate(diff_energy);

        if (m_calculatePressure)
        {
            const double newPressureParticle {pressureParticle(m_temp, indexTranslation, positionParticleTranslation, m_positionArray, neighborIList, m_diameterArray, m_squareRc, m_lengthCube, m_halfLengthCube)};
            const double oldPressureParticle {pressureParticle(m_temp, indexTranslation, m_positionArray[indexTranslation], m_positionArray, neighborIList, m_diameterArray, m_squareRc, m_lengthCube, m_halfLengthCube)};
            m_pressure += newPressureParticle - oldPressureParticle;
        }

        m_stepDisplacementMatrix[indexTranslation] = vectorSum(m_stepDisplacementMatrix[indexTranslation],
                                                               randomVector);
        m_positionArray[indexTranslation] = positionParticleTranslation;
    }
}

/*******************************************************************************
 * This function returns a tentative new particle position.
 *
 * @param indexTranslation Translated particle's index.
 *        randomVector Random vector of size 3. The components of the vector are
 *                     taken from a uniform distribution U(-m_rBox, m_rBox).
 * @return Tentative new particle position.
 ******************************************************************************/
void MonteCarlo::mcSwap()
{
    int indexSwap1 { randomIntGenerator(0, m_nParticles - 1) };
    int indexSwap2;
    // int indexSwap2 { randomIntGenerator(0, m_nParticles - 1) };
    // !!!! THIS NEXT PART IS ONLY BECAUSE WE STUDY TRI-MERS VERY SPECIFIC!!!
    indexSwap1 -= indexSwap1 % 3;
    indexSwap2 = indexSwap1 + 2;

    const double sigma1 { m_diameterArray[indexSwap1] };
    const double sigma2 { m_diameterArray[indexSwap2] };
    double energyParticle1;
    double energyParticle2;
    double energyParticleSwap1;
    double energyParticleSwap2;
    const std::vector<int> neighborIList1 { m_neighborList[indexSwap1] };
    const std::vector<int> neighborIList2 { m_neighborList[indexSwap2] };


    if (m_simulationMol == "polymer")
    {

        const std::vector<int> bondsI1 { m_bondsMatrix[indexSwap1] };
        const std::vector<int> bondsI2 { m_bondsMatrix[indexSwap2] };

        energyParticle1 = energyParticlePolymer(indexSwap1, m_positionArray[indexSwap1],
                                                m_positionArray, neighborIList1, m_diameterArray, bondsI1,
                                                m_squareRc, m_lengthCube, m_halfLengthCube, m_squareR0, m_feneK, m_bondType);

        energyParticle2 = energyParticlePolymer(indexSwap2, m_positionArray[indexSwap2],
                                                m_positionArray, neighborIList2, m_diameterArray, bondsI2,
                                                m_squareRc, m_lengthCube, m_halfLengthCube, m_squareR0, m_feneK, m_bondType, indexSwap1);

        m_diameterArray[indexSwap1] = sigma2;
        m_diameterArray[indexSwap2] = sigma1;

        energyParticleSwap1 = energyParticlePolymer(indexSwap1, m_positionArray[indexSwap1],
                                                    m_positionArray, neighborIList1, m_diameterArray, bondsI1,
                                                    m_squareRc, m_lengthCube, m_halfLengthCube, m_squareR0, m_feneK, m_bondType);

        energyParticleSwap2 = energyParticlePolymer(indexSwap2, m_positionArray[indexSwap2],
                                                    m_positionArray, neighborIList2, m_diameterArray, bondsI2,
                                                    m_squareRc, m_lengthCube, m_halfLengthCube, m_squareR0, m_feneK, m_bondType, indexSwap1);

    }
    else
    {
        energyParticle1 = energyParticle(indexSwap1, m_positionArray[indexSwap1],
                                         m_positionArray, neighborIList1, m_diameterArray,
                                         m_squareRc, m_lengthCube, m_halfLengthCube);

        energyParticle2 = energyParticle(indexSwap2, m_positionArray[indexSwap2],
                                         m_positionArray, neighborIList2, m_diameterArray,
                                         m_squareRc, m_lengthCube, m_halfLengthCube, indexSwap1);

        m_diameterArray[indexSwap1] = sigma2;
        m_diameterArray[indexSwap2] = sigma1;

        energyParticleSwap1 = energyParticle(indexSwap1, m_positionArray[indexSwap1],
                                             m_positionArray, neighborIList1, m_diameterArray,
                                             m_squareRc, m_lengthCube, m_halfLengthCube);

        energyParticleSwap2 = energyParticle(indexSwap2, m_positionArray[indexSwap2],
                                             m_positionArray, neighborIList2, m_diameterArray,
                                             m_squareRc, m_lengthCube, m_halfLengthCube,indexSwap1);
    }
    const double diff_energy{ energyParticleSwap1 + energyParticleSwap2 - energyParticle1 - energyParticle2 };
    // Metropolis criterion
    const bool acceptMove { metropolis( diff_energy) };

    // If the move is accepted, then the energy, the position array and the displacement array can be updated.
    // If the m_calculatePressure is set to True, then the pressure is calculated.
    if ( acceptMove )
    {
        generalUpdate( diff_energy );
        m_acceptanceRateSwap += 1. / m_nParticles;
        if (m_calculatePressure)
        {
            const double pressureParticleSwap1 {pressureParticle(m_temp, indexSwap1,
                                                           m_positionArray[indexSwap1],
                                                           m_positionArray, neighborIList1,
                                                           m_diameterArray, m_squareRc,
                                                           m_lengthCube, m_halfLengthCube)};

            const double pressureParticleSwap2 {pressureParticle(m_temp, indexSwap2,
                                                           m_positionArray[indexSwap2],
                                                           m_positionArray, neighborIList2,
                                                           m_diameterArray, m_squareRc,
                                                           m_lengthCube, m_halfLengthCube)};

            m_diameterArray[indexSwap1] = sigma1;
            m_diameterArray[indexSwap2] = sigma2;

            const double pressureParticle1 {pressureParticle(m_temp, indexSwap1,
                                                           m_positionArray[indexSwap1],
                                                           m_positionArray, neighborIList1,
                                                           m_diameterArray, m_squareRc,
                                                           m_lengthCube, m_halfLengthCube)};


            const double pressureParticle2 {pressureParticle(m_temp, indexSwap2,
                                                           m_positionArray[indexSwap2],
                                                           m_positionArray, neighborIList2,
                                                           m_diameterArray, m_squareRc,
                                                           m_lengthCube, m_halfLengthCube)};

            m_diameterArray[indexSwap1] = sigma2;
            m_diameterArray[indexSwap2] = sigma1;

            m_pressure += pressureParticleSwap1 + pressureParticleSwap2 - pressureParticle1 - pressureParticle2;
        }
    }
    else
    {
        m_diameterArray[indexSwap1] = sigma1;
        m_diameterArray[indexSwap2] = sigma2;
    }
}


void MonteCarlo::generalUpdate(double diff_energy)
{

    const double newEnergy {m_energy + diff_energy};
    m_energy = newEnergy;
    m_acceptanceRate += 1. / m_nParticles; // increment of the acceptance rate.

}


std::vector<int> MonteCarlo::findNeighborIList(int indexTranslation)
{
    std::vector<int> neighborIList { };

    if ( m_neighMethod == "verlet" )
    {
        neighborIList = m_neighborList[indexTranslation]; // particle's neighbor list
    }

    else
    {
        neighborIList.resize ( m_nParticles );
        std::iota (std::begin(neighborIList), std::end(neighborIList), 0);
    }
    return neighborIList;
}

/*******************************************************************************
 * This function is an implementation of the Metropolis algorithm.
 * Two energies are compared: the energy of the new configuration (after a
 * Monte Carlo move) and the energy of the former configuration.
 * The Metropolis algorithm decides if the MC move is accepted or not according
 * to the following conditions:
 * - if newEnergy < energy then the move is accepted.
 * - if newEnergy > energy the the move is accepted at a certain probability
 *   proportional to exp((energy-newEnergy)/kT).
 *
 * @param diff_energy New tentative configuration's energy minus old configuration's energy.
 *
 * @return Returns true if the move is accepted and False otherwise.
 ******************************************************************************/
bool MonteCarlo::metropolis(double diff_energy) const
{
	if (diff_energy < 0)
		return true;

	else
	{
		const double randomDouble { randomDoubleGenerator(0., 1.) } ;
		const double threshold { exp(( - diff_energy ) / m_temp) }; // we consider k=1
		return threshold > randomDouble;
	}
}

/*******************************************************************************
 * This function checks if the maximum displacement of one particle is superior
 * than m_squareRDiff. If this is true, then the neighbor list is updated and
 * m_interDisplacementMatrix is reinitialized to zero. The criterion shouldn't
 * allow for any errors?
 ******************************************************************************/
void MonteCarlo::checkStepDisplacement()
{
   const  std::vector<double> squareDispVector = getSquareNormRowMatrix(m_interDisplacementMatrix);

    for (int i=0; i < m_nParticles; i++)
    {
        const double max_rc {( m_maxDiam + m_diameterArray[i]) / 2 * m_rC};
        const double thresh = { pow ( (m_rSkin - max_rc) / 2, 2)};
        if ( squareDispVector[i] > thresh)
        {
            createNeighborList();
            std::fill(m_interDisplacementMatrix.begin(), m_interDisplacementMatrix.end(), std::vector<double>(3, 0));
            break;
        }
    }
}