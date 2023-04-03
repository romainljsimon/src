/*
 * MonteCarlo.cpp
 *
 *  Created on: 27 oct. 2022
 *      Author: rsimon
 */

#include <iostream>
#include <utility>
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


MonteCarlo::MonteCarlo ( std::string simulationMol, std::vector<std::vector<double>> positionArray,
						 std::vector<double> diameterArray, std::vector<int> moleculeType,const double rc,
						 const double lengthCube, const double temp,
						 const double rbox, const double rskin, const int saveUpdate,
						 std::string  folderPath, std::string  neighMethod,
						 const int timeSteps, const double r0, const double feneK)

	: m_simulationMol (std::move(simulationMol ))
	, m_positionArray (std::move( positionArray ))
	, m_diameterArray (std::move( diameterArray ))
	, m_moleculeType (std::move( moleculeType ))
	, m_squareRc { pow ( rc, 2 ) }
	, m_lengthCube { lengthCube }
	, m_temp { temp }
	, m_rbox { rbox }
	, m_squareRskin { pow ( rskin, 2 ) }
	, m_saveUpdate { saveUpdate }
	, m_folderPath (std::move( folderPath ))
	, m_neighMethod (std::move( neighMethod ))
	, m_timeSteps { timeSteps }
	, m_squareR0 { pow ( r0, 2 ) }
	, m_feneK { feneK }
	, m_squareRdiff{ pow ( (rskin - rc) / 2, 2)}

{
	m_nParticles = static_cast<int>( m_positionArray.size() );
	m_interDisplacementMatrix.resize(m_nParticles, std::vector<double>(3, 0));
	m_totalDisplacementMatrix.resize(m_nParticles, std::vector<double>(3, 0));
	m_stepDisplacementMatrix.resize(m_nParticles, std::vector<double>(3, 0));
	m_neighborList.resize(m_nParticles);

	createNeighborList();

	if (m_simulationMol == "polymer")
	{
		m_bondsMatrix = readBondsTXT(m_folderPath + "/bonds.txt");
		m_energy = energySystemPolymer(m_positionArray, m_diameterArray, m_bondsMatrix,
									   m_neighborList, m_squareRc, m_lengthCube, m_squareR0, m_feneK);
		// m_energy = 0;
	}

	else
	{
		// double density {m_nParticles / pow(m_lengthCube, 3.)};
		m_energy = energySystem(m_positionArray, m_diameterArray, m_neighborList, m_squareRc, m_lengthCube);  // + 8. / 3 * density * 3.14 * m_nParticles * 1 / pow(m_squareRc, 3./2) * (1. / 3 * pow(1 / m_squareRc, 3. ) - 1);
		// double corrPressure {16. * 3.14 / 3. * density / ( m_temp *  pow(m_squareRc, 3./2)) * (2. / (3 * pow(m_squareRc, 3. )) - 1.) };

		if (m_calculatePressure)
		{
			m_pressure = pressureSystem(m_temp, m_positionArray, m_diameterArray, m_squareRc, m_lengthCube); //+ corrPressure;
		}
	}
}


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
	std::string prename (m_folderPath + "/outXYZ/position");
	std::string extname {".xyz"};
	std::string prenameDisp (m_folderPath + "/disp/displacement");
	std::string extnameDisp {".txt"};
    std::vector<double> radiusArray (divideVectorByScalar(m_diameterArray, 2));
	saveInXYZ(m_positionArray,  radiusArray, m_moleculeType, m_lengthCube, prename + std::to_string(0) + extname );
	saveDoubleTXT(m_energy / m_nParticles, m_folderPath + "/outE.txt");
	saveDisplacement(m_totalDisplacementMatrix, prenameDisp + std::to_string(0) + extnameDisp);

	std::vector<int> saveTimeStepArray ( createSaveTime(m_timeSteps, m_saveUpdate, 1.5));

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
            saveInXYZ ( m_positionArray,  radiusArray, m_moleculeType, m_lengthCube, prename + std::to_string(i+1) + extname );
			saveDisplacement ( m_totalDisplacementMatrix, prenameDisp + std::to_string(i+1) + extnameDisp );
			++save_index;
		}

		if (i % 50 == 0)
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
	saveInXYZ(m_positionArray,  radiusArray, m_moleculeType, m_lengthCube, prename + std::to_string(m_timeSteps) + extname );
	saveDisplacement(m_totalDisplacementMatrix, prenameDisp + std::to_string(m_timeSteps) + extnameDisp);
	saveDoubleTXT( m_errors, m_folderPath + "/errors.txt");

	std::cout << "MC move acceptance rate: " << m_acceptanceRate << "\n";
	std::cout << "Neighbor list update rate: " << m_updateRate / m_timeSteps << "\n";
	std::cout << "Number of neighbor list errors: " << m_errors << "\n";

}


/*******************************************************************************
 * Creation or update of the Verlet neighbor list. The neighbor list is
 * implemented as an 2d array suchs as row i of the array is particle i's
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
	++m_updateRate;
	std::vector<std::vector<int>> oldNeighborList = m_neighborList;
	m_neighborList.clear();
	m_neighborList.resize(m_nParticles);

	for (int i = 0; i < m_nParticles - 1; i++)
	{
		std::vector<double> positionParticle = m_positionArray[i];

		for (int j = i + 1; j < m_nParticles; j++)
		{

			double squareDistance { squareDistancePair (positionParticle, m_positionArray[j], m_lengthCube) };

			if (squareDistance < m_squareRskin)
			{
				// j and i are neighbors

				m_neighborList[i].push_back(j); // j is on row i
				m_neighborList[j].push_back(i); // i is on row j

				if (static_cast<int>(oldNeighborList[i].size()) > 0) // If this size is 0, then it is the first time the neighbor list is calculated.
				{

					// Compare the new neighbor list with the new neighbor list

					if (!std::binary_search(oldNeighborList[i].begin(), oldNeighborList[i].end(), j))
					{
						if (squareDistance < m_squareRc)
						{
							++m_errors;
						}
					}
				}

			}
		}

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
	int indexTranslation { randomIntGenerator(0, m_nParticles - 1) }; // randomly chosen particle
	std::vector<double> randomVector ( randomVectorDoubleGenerator(3, -m_rbox, m_rbox) );
	std::vector<double> positionParticleTranslation = mcTranslation ( indexTranslation, randomVector );
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

	double oldEnergyParticle {};
	double newEnergyParticle {};

	if (m_simulationMol == "polymer")
	{
 		std::vector<int> bondsI ( m_bondsMatrix[indexTranslation] );

		oldEnergyParticle = energyParticlePolymer (indexTranslation, m_positionArray[indexTranslation],
													m_positionArray, neighborIList, m_diameterArray, bondsI,
													m_squareRc, m_lengthCube, m_squareR0, m_feneK);
		newEnergyParticle = energyParticlePolymer (indexTranslation, positionParticleTranslation,
													m_positionArray, neighborIList, m_diameterArray, bondsI,
													m_squareRc, m_lengthCube, m_squareR0, m_feneK);

	}

	else
	{
		oldEnergyParticle = energyParticle (indexTranslation, m_positionArray[indexTranslation], m_positionArray,
											neighborIList, m_diameterArray, m_squareRc, m_lengthCube);
		newEnergyParticle = energyParticle (indexTranslation, positionParticleTranslation, m_positionArray,
											neighborIList, m_diameterArray, m_squareRc, m_lengthCube);
	}

	// Metropolis criterion
	bool acceptMove { metropolis (newEnergyParticle, oldEnergyParticle) } ;

	// If the move is accepted, then the energy, the position array and the displacement array can be updated.
	// If the m_calculatePressure is set to True, then the pressure is calculated.
	if ( acceptMove )
	{

		if (m_calculatePressure)
		{
			double newPressureParticle {pressureParticle(m_temp, indexTranslation, positionParticleTranslation, m_positionArray, neighborIList, m_diameterArray, m_squareRc, m_lengthCube)};
			double oldPressureParticle {pressureParticle(m_temp, indexTranslation, m_positionArray[indexTranslation], m_positionArray, neighborIList, m_diameterArray, m_squareRc, m_lengthCube)};
			m_pressure += newPressureParticle - oldPressureParticle;
		}

		double newEnergy {m_energy - oldEnergyParticle + newEnergyParticle};
		m_energy = newEnergy;
		m_stepDisplacementMatrix[indexTranslation] = vectorSum(m_stepDisplacementMatrix[indexTranslation], randomVector);
		m_positionArray[indexTranslation] = positionParticleTranslation;
		m_acceptanceRate += 1. / m_nParticles; // increment of the acceptance rate.

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
	std::vector<double> randomVector ( randomVectorDoubleGenerator(3, -m_rbox, m_rbox) );

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
 *                     taken from a uniform distribution U(-m_rbox, m_rbox).
 * @return Tentative new particle position.
 ******************************************************************************/
std::vector<double> MonteCarlo::mcTranslation(int indexTranslation, const std::vector<double>& randomVector)
{

	std::vector<double> positionParticleTranslation = m_positionArray[indexTranslation];
	positionParticleTranslation = vectorSum (positionParticleTranslation, randomVector);
	positionParticleTranslation = periodicBC (positionParticleTranslation, m_lengthCube);
	return positionParticleTranslation;
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
 * @param newEnergy New tentative configuration's energy.
 *        energy Old configuration's energy.
 *
 * @return Returns true if the move is accepted and False otherwise.
 ******************************************************************************/
bool MonteCarlo::metropolis(double newEnergy, double energy)
{
	if (newEnergy < energy)
		return true;

	else
	{
		double randomDouble { randomDoubleGenerator(0., 1.) } ;
		double threshold { exp((energy - newEnergy) / m_temp) }; // we consider k=1
		return threshold > randomDouble;
	}
}

/*******************************************************************************
 * This function checks if the maximum displacement of one particle is superior
 * than m_squareRdiff. If this is true, then the neighbor list is updated and
 * m_interDisplacementMatrix is reinitialized to zero. The criterion shouldn't
 * allow for any errors?
 ******************************************************************************/
void MonteCarlo::checkStepDisplacement()
{
	std::vector<double> squareDispVector = getSquareNormRowMatrix(m_interDisplacementMatrix);

	if (getMaxVector ( squareDispVector ) > m_squareRdiff)
	{
		createNeighborList();
		std::fill(m_interDisplacementMatrix.begin(), m_interDisplacementMatrix.end(), std::vector<double>(3, 0));
	}
}
