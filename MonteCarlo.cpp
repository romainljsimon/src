/*
 * MonteCarlo.cpp
 *
 *  Created on: 27 oct. 2022
 *      Author: rsimon
 */

#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>

#include <string>
#include <numeric>
#include <algorithm>
#include "MonteCarlo.h"
#include "random.h"
#include "energy.h"
#include "readSaveFile.h"
#include "util.h"



MonteCarlo::MonteCarlo ( std::string simulationMol, std::vector<std::vector<double>> positionArray,
						 std::vector<double> radiusArray, std::vector<int> moleculeType,const double rc,
						 const double lengthCube, const double temp,
						 const double rbox, const double rskin, const int neighUpdate,
						 const std::string folderPath, const std::string neighMethod,
						 const int timeSteps, const double r0, const double feneK)

	: m_simulationMol ( simulationMol )
	, m_positionArray ( positionArray )
	, m_radiusArray ( radiusArray )
	, m_moleculeType ( moleculeType )
	, m_squareRc { pow ( rc, 2 ) }
	, m_lengthCube { lengthCube }
	, m_temp { temp }
	, m_rbox { rbox }
	, m_squareRskin { pow ( rskin, 2 ) }
	, m_neighUpdate { neighUpdate }
	, m_folderPath ( folderPath )
	, m_neighMethod ( neighMethod )
	, m_timeSteps { timeSteps }
	, m_squareR0 { pow ( r0, 2 ) }
	, m_feneK { feneK }
	, m_squareRdiff{ pow ( (rskin - rc) / 2, 2)}

{
	m_nParticles = static_cast<int>( m_positionArray.size() );
	m_interDisplacementMatrix.resize(m_nParticles, std::vector<double>(3, 0));
	m_totalDisplacementMatrix.resize(m_nParticles, std::vector<double>(3, 0));
	m_stepDisplacementMatrix.resize(m_nParticles, std::vector<double>(3, 0));

	//		m_cellList.resize(m_nParticles,
	//		std::vector<std::vector<std::vector<int>>>(m_nParticles,
	//		std::vector<std::vector<int>>(m_nParticles)));
	m_neighborList.resize(m_nParticles);

	if (m_simulationMol == "polymer")
	{
		m_bondsMatrix = readBondsTXT(m_folderPath + "/bonds.txt");
		m_energy = energySystemPolymer(m_positionArray, m_radiusArray, m_bondsMatrix,
									   m_squareRc, m_lengthCube, m_squareR0, m_feneK);
	}

	else
	{
		// double density {m_nParticles / pow(m_lengthCube, 3.)};
		m_energy = energySystem(m_positionArray, m_radiusArray, m_squareRc, m_lengthCube);  // + 8. / 3 * density * 3.14 * m_nParticles * 1 / pow(m_squareRc, 3./2) * (1. / 3 * pow(1 / m_squareRc, 3. ) - 1);
		// double corrPressure {16. * 3.14 / 3. * density / ( m_temp *  pow(m_squareRc, 3./2)) * (2. / (3 * pow(m_squareRc, 3. )) - 1.) };
		if (m_calculatePressure)
		{
			m_pressure = pressureSystem(m_temp, m_positionArray, m_radiusArray, m_squareRc, m_lengthCube); //+ corrPressure;
		}
	}
	createNeighborList();
}



void MonteCarlo::mcTotal()

/* This function is the core of  Monte Carlo program. It iterates over a certain amount of steps that the user can choose.
 * At each step a Monte Carlo move is made. For now the Monte Carlo move is a translation.
 * The position of each particle is written in a .xyz file at the end of each step.
 */

{
	std::string prename (m_folderPath + "/outXYZ/position");
	std::string extname {".xyz"};
	std::string prenameDisp (m_folderPath + "/disp/displacement");
	std::string extnameDisp {".txt"};
	saveInXYZ(m_positionArray,  m_radiusArray, m_moleculeType, m_lengthCube, prename + std::to_string(0) + extname );
	saveEnergyTXT(m_energy / m_nParticles, m_folderPath + "/outE.txt");
	int thresh {1};
	for (int i = 5000000; i < 5000000 + m_timeSteps; i++)
	{

		if (m_calculatePressure)
		{
			std::ofstream outPressure;
			outPressure.open(m_folderPath + "/pressure.txt", std::ios_base::app);
			outPressure << m_pressure;
			outPressure << "\n";
			outPressure.close();
		}

		for (int j = 0; j < m_nParticles; j++)
		{
			mcMove();
		}
		m_interDisplacementMatrix = matrixSum(m_interDisplacementMatrix, m_stepDisplacementMatrix);
		m_totalDisplacementMatrix = matrixSum(m_totalDisplacementMatrix, m_stepDisplacementMatrix);
		std::fill(m_stepDisplacementMatrix.begin(), m_stepDisplacementMatrix.end(), std::vector<double>(3, 0));

		if (m_neighMethod == "verlet")
		{
			checkStepDisplacement();
		}
		if ((i % 1000) == 0)
		{
			saveInXYZ(m_positionArray,  m_radiusArray, m_moleculeType, m_lengthCube, prename + std::to_string(i) + extname );
		}

		if ((i % thresh) == 0)
		{
			//saveInXYZ(m_positionArray,  m_radiusArray, m_moleculeType, m_lengthCube, prename + std::to_string(i) + extname );
			saveDisplacement(m_totalDisplacementMatrix, prenameDisp + std::to_string(i) + extnameDisp);
		}
		if (i > (10 * thresh - 1))
		{
			thresh = 10 * thresh;
		}
		saveEnergyTXT(m_energy / m_nParticles, m_folderPath + "/outE.txt");


	}

	std::ofstream outErrors;
	outErrors.open(m_folderPath + "/errors.txt", std::ios_base::app);
	outErrors << m_errors;
	outErrors << "\n";
	outErrors.close();
	std::cout << "The acceptance rate is: " << m_acceptanceRate / ( m_nParticles * m_timeSteps ) << "\n";
	std::cout << "The neighbor list update rate is: " << m_updateRate / m_timeSteps << "\n";

}


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
				if (static_cast<int>(oldNeighborList[i].size()) > 0)
				{
					if (!std::binary_search(oldNeighborList[i].begin(), oldNeighborList[i].end(), j))
					{
						if (squareDistance < m_squareRc)
						{
							++m_errors;
						}
					}
				}
				m_neighborList[i].push_back(j);
				m_neighborList[j].push_back(i);

			}
		}

	}
}

//void MonteCarlo::createCellAndIndexList()
//{
//
//	// int cellListSize {static_cast<int>(m_lengthCube / m_rskin)};
//
//	for (int i = 0; i < m_nParticles - 1; i++)
//	{
//		std::vector<double> positionParticle = m_positionArray[i];
//		std::vector<double> positionBoxLength = divideVectorByScalar(positionParticle, m_rskin);
//		std::vector<int> boxIndex ( positionBoxLength.begin(), positionBoxLength.end() );
//		m_cellList[boxIndex[0]][boxIndex[1]][boxIndex[2]].push_back(i);
//		m_particleIndexCell[i] = boxIndex;
//	}
//
//}

void MonteCarlo::mcMove()
{
	/*
	 *This function implements a Monte Carlo move: translation of a random particle, calculation of the energy of the new system and then acceptation or not of the move.
	 */
	int indexTranslation { randomIntGenerator(0, m_nParticles - 1) };
	std::vector<double> randomVector ( randomVectorDoubleGenerator(3, -m_rbox, m_rbox) );

	std::vector<double> positionParticleTranslation = mcTranslation ( indexTranslation, randomVector );
	std::vector<int> neighborIList { };

	if (m_neighMethod == "verlet")
	{
		neighborIList = m_neighborList[indexTranslation];
	}

	else
	{
		neighborIList.resize (m_nParticles);
		std::iota (std::begin(neighborIList), std::end(neighborIList), 0);
	}

	double oldEnergyParticle {};
	double newEnergyParticle {};

	if (m_simulationMol == "polymer")
	{
 		std::vector<int> bondsI ( m_bondsMatrix[indexTranslation] );

		oldEnergyParticle = energyParticlePolymer (indexTranslation, m_positionArray[indexTranslation],
													m_positionArray, neighborIList, m_radiusArray, bondsI,
													m_squareRc, m_lengthCube, m_squareR0, m_feneK);
		newEnergyParticle = energyParticlePolymer (indexTranslation, positionParticleTranslation,
													m_positionArray, neighborIList, m_radiusArray, bondsI,
													m_squareRc, m_lengthCube, m_squareR0, m_feneK);

	}

	else
	{
		oldEnergyParticle = energyParticle (indexTranslation, m_positionArray[indexTranslation], m_positionArray,
											neighborIList, m_radiusArray, m_squareRc, m_lengthCube);
		newEnergyParticle = energyParticle (indexTranslation, positionParticleTranslation, m_positionArray,
											neighborIList, m_radiusArray, m_squareRc, m_lengthCube);
	}

	bool acceptMove { metropolis (newEnergyParticle, oldEnergyParticle) } ;

//	if ((indexTranslation==0) || (indexTranslation==(m_nParticles-1) || (indexTranslation==40) || (indexTranslation==39) ))
//	{
//		acceptMove=false;
//	}

	if (acceptMove)
	{
		if (m_calculatePressure)
		{
			double newPressureParticle {pressureParticle(m_temp, indexTranslation, positionParticleTranslation, m_positionArray, neighborIList, m_radiusArray, m_squareRc, m_lengthCube)};
			double oldPressureParticle {pressureParticle(m_temp, indexTranslation, m_positionArray[indexTranslation], m_positionArray, neighborIList, m_radiusArray, m_squareRc, m_lengthCube)};
			m_pressure += newPressureParticle - oldPressureParticle;
		}

		double newEnergy {m_energy - oldEnergyParticle + newEnergyParticle};
		m_energy = newEnergy;
		m_stepDisplacementMatrix[indexTranslation] = vectorSum(m_stepDisplacementMatrix[indexTranslation], randomVector);
		m_positionArray[indexTranslation] = positionParticleTranslation;

		++m_acceptanceRate;

	}

}

std::vector<double> MonteCarlo::mcTranslation(int indexTranslation, std::vector<double> randomVector)
{
	/*
	 *This function chooses a particle randomly and translates it randomly inside of box of size rbox.
	 */

	std::vector<double> positionParticleTranslation = m_positionArray[indexTranslation];
	positionParticleTranslation = vectorSum (positionParticleTranslation, randomVector);
	positionParticleTranslation = periodicBC (positionParticleTranslation, m_lengthCube);
	return positionParticleTranslation;
}


bool MonteCarlo::metropolis(double newEnergy, double energy)
/*
 *This function is an implementation of the Metropolis algorithm.
 *Two energies are compared: the energy of the new configuration (after a mcMove) and the energy of former configuration.
 *The Metropolis algorithm decides if the MC move is accepted or not according to the following conditions:
 *- if newEnergy < energy then the move is accepted
 *- if newEnergy > energy the the move is accepted at a certain probability proportional to exp((energy-newEnergy)/kT)
 */
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

void MonteCarlo::checkStepDisplacement()
{
	std::vector<double> squareDispVector = getSquareNormRowMatrix(m_interDisplacementMatrix);

	if (getMaxVector ( squareDispVector ) > m_squareRdiff)
	{
		createNeighborList();
		std::fill(m_interDisplacementMatrix.begin(), m_interDisplacementMatrix.end(), std::vector<double>(3, 0));
	}
}
