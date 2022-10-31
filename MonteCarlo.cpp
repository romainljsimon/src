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



MonteCarlo::MonteCarlo ( std::vector<std::vector<double>> positionArray, std::vector<double> radiusArray,
						double rc, double lengthCube, double temp, double rbox, double rskin, int neighUpdate,
						std::string folderPath, std::string neighMethod, int numberIteration )
    : m_positionArray { positionArray }
	, m_radiusArray {radiusArray}
	, m_rc {rc}
	, m_lengthCube { lengthCube }
	, m_temp { temp }
	, m_rbox { rbox }
	, m_rskin { rskin }
	, m_neighUpdate { neighUpdate }
	, m_folderPath { folderPath }
	, m_neighMethod { neighMethod }
	, m_numberIteration { numberIteration }
{
	m_nParticles = static_cast<int>( m_positionArray.size() );
	m_particleIndexCell.resize(m_nParticles);
	m_cellList.resize(m_nParticles,
			std::vector<std::vector<std::vector<int>>>(m_nParticles,
			std::vector<std::vector<int>>(m_nParticles)));
	m_neighborList.resize(m_nParticles);

}



void MonteCarlo::mcTotal()

/* This function is the core of  Monte Carlo program. It iterates over a certain amount of steps that the user can choose.
 * At each step a Monte Carlo move is made. For now the Monte Carlo move is a translation.
 * The position of each particle is written in a .xyz file at the end of each step.
  */

{
	m_energy = energySystem(m_positionArray, m_radiusArray, m_rc, m_lengthCube);
	createNeighborList();
	createCellAndIndexList();
	std::string prename {m_folderPath + "/outXYZ/position"};
	std::string extname {".xyz"};

	for (int i = 1; i < m_numberIteration + 1; i++)
	{
		mcMove();
		if (((i % m_neighUpdate) == 0) && (m_neighMethod == "verlet") )
		{
			createNeighborList();
		}
		if (i > m_numberIteration - 10)
		{
			saveInXYZ(m_positionArray,  m_radiusArray, prename + std::to_string(i) + extname );
		}

		saveEnergyTXT(m_energy / m_nParticles, m_folderPath + "/outE.txt");

	}
	std::ofstream outErrors;
	outErrors.open(m_folderPath + "/errors.txt", std::ios_base::app);
	outErrors << m_errors;
	outErrors << "\n";
	outErrors.close();
}

//void MonteCarlo::createNeighborList()
//{
//	m_neighborList.clear();
//	m_neighborList.resize(m_nParticles);
//
//	double squareSkin {pow(m_rskin, 2)};
//
//	for (int i = 0; i < m_nParticles - 1; i++)
//	{
//		std::vector<double> positionParticle = m_positionArray[i];
//
//		for (int j = i + 1; j < m_nParticles; j++)
//		{
//
//			if (squareDistancePair (positionParticle, m_positionArray[j], m_lengthCube) < squareSkin)
//			{
//
//				m_neighborList[i].push_back(j);
//				m_neighborList[j].push_back(i);
//			}
//		}
//
//	}
//}

void MonteCarlo::createNeighborList()
{
	std::vector<std::vector<int>> oldNeighborList = m_neighborList;
	m_neighborList.clear();
	m_neighborList.resize(m_nParticles);

	double squareSkin {pow(m_rskin, 2)};

	for (int i = 0; i < m_nParticles - 1; i++)
	{
		std::vector<double> positionParticle = m_positionArray[i];

		for (int j = i + 1; j < m_nParticles; j++)
		{

			double squareDistance { squareDistancePair (positionParticle, m_positionArray[j], m_lengthCube) };

			if (squareDistance < squareSkin)
			{
				if (static_cast<int>(oldNeighborList[i].size()) > 0)
				{
					if(std::find(oldNeighborList[i].begin(), oldNeighborList[i].end(), j) != oldNeighborList[i].end())
					{
						;
					}
					else
					{
						if (squareDistance < pow(m_rc, 2))
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

void MonteCarlo::createCellAndIndexList()
{

	int cellListSize {static_cast<int>(m_lengthCube / m_rskin)};

	for (int i = 0; i < m_nParticles - 1; i++)
	{
		std::vector<double> positionParticle = m_positionArray[i];
		std::vector<double> positionBoxLength = divideVectorByScalar(positionParticle, m_rskin);
		std::vector<int> boxIndex ( positionBoxLength.begin(), positionBoxLength.end() );
		m_cellList[boxIndex[0]][boxIndex[1]][boxIndex[2]].push_back(i);
		m_particleIndexCell[i] = boxIndex;
	}

}

void MonteCarlo::mcMove()
{
	/*
	 *This function implements a Monte Carlo move: translation of a random particle, calculation of the energy of the new system and then acceptation or not of the move.
	 */
	mcTranslation();
	std::vector<int> neighborIList { };
	if (m_neighMethod == "verlet")
	{
		neighborIList = m_neighborList[m_indexTranslation];
	}
	else
	{
		neighborIList.resize (m_nParticles);
		std::iota (std::begin(neighborIList), std::end(neighborIList), 0);
	}
	double oldEnergyParticle { energyParticle(1, m_indexTranslation, m_positionArray[m_indexTranslation], m_positionArray, neighborIList, m_radiusArray, m_rc, m_rskin, m_lengthCube) };
	double newEnergyParticle { energyParticle(0, m_indexTranslation, m_positionParticleTranslation, m_positionArray, neighborIList, m_radiusArray, m_rc, m_rskin, m_lengthCube) };
	metropolis (newEnergyParticle, oldEnergyParticle);

	if (m_acceptMove)
	{
		double newEnergy {m_energy - oldEnergyParticle + newEnergyParticle};
		m_energy = newEnergy;
		m_positionArray[m_indexTranslation] = m_positionParticleTranslation;
	}
	m_positionParticleTranslation.clear();

}

void MonteCarlo::mcTranslation()
{
	/*
	 *This function chooses a particle randomly and translates it randomly inside of box of size rbox.
	 */

	m_indexTranslation =  randomIntGenerator(0, m_nParticles-1);
	m_positionParticleTranslation = m_positionArray[m_indexTranslation];
	std::vector<double> randomVector ( randomVectorDoubleGenerator(3, -m_rbox, m_rbox) );
	m_positionParticleTranslation = vectorSum (m_positionParticleTranslation, randomVector);
	m_positionParticleTranslation = periodicBC (m_positionParticleTranslation, m_lengthCube);

}


void MonteCarlo::metropolis(double newEnergy, double energy)
/*
 *This function is an implementation of the Metropolis algorithm.
 *Two energies are compared: the energy of the new configuration (after a mcMove) and the energy of former configuration.
 *The Metropolis algorithm decides if the MC move is accepted or not according to the following conditions:
 *- If newEnergy < energy then the move is accepted
 *- If newEnergy > energy the the move is accepted at a certain probability proportional to exp((energy-newEnergy)/kT)
 */
{
	if (newEnergy < energy)
		m_acceptMove = true;

	else
	{
		double randomDouble { randomDoubleGenerator(0., 1.) } ;
		double threshold { exp((energy - newEnergy) / m_temp) }; // we consider k=1
		m_acceptMove = threshold > randomDouble;
	}
}

