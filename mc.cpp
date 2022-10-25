/*
 * mc.cpp
 *
 *  Created on: 4 oct. 2022
 *      Author: romainsimon
 */

#include <iostream>
#include <vector>
#include <math.h>
#include <string>

#include "random.h"
#include "energy.h"
#include "readSaveFile.h"
#include "util.h"

struct enePos
{
	double ene;
	std::vector<std::vector<double>> posMatrix;
};

struct indexPosParticle
{
	int index;
	std::vector<double> posParticle;
};

bool metropolis(double newEnergy, double energy, double temp)
/*
 *This function is an implementation of the Metropolis algorithm.
 *Two energies are compared: the energy of the new configuration (after a mcMove) and the energy of former configuration.
 *The Metropolis algorithm decides if the MC move is accepted or not according to the following conditions:
 *- If newEnergy < energy then the move is accepted
 *- If newEnergy > energy the the move is accepted at a certain probability proportional to exp((energy-newEnergy)/kT)
 */
{
	if (newEnergy < energy)
		return true;

	else
	{
		double randomDouble { randomDoubleGenerator(0., 1.) } ;
		double threshold { exp((energy - newEnergy) / temp) }; // we consider k=1
		return threshold > randomDouble;
	}
}

std::vector<double> periodicBC(std::vector<double> positionParticle, double lengthCube)
/*
 *This function is an implementation of the periodic Boundary conditions.
 *If a particle gets out of the simulation box from one of the sides, it gets back in the box from the opposite side.
 */
{
	int positionParticleSize { static_cast<int>(positionParticle.size()) };
	for (int i = 0; i < positionParticleSize; i++)
	{
		if (positionParticle[i]< 0)
			positionParticle[i] += lengthCube;

		else if (positionParticle[i] > lengthCube)
			positionParticle[i] -= lengthCube;
	}
	return positionParticle;
}

indexPosParticle mcTranslation(std::vector<std::vector<double>>& positionArray, double rbox, double lengthCube)
{
	/*
	 *This function chooses a particle randomly and translates it randomly inside of box of size rbox.
	 */
	indexPosParticle translation;

	int i { randomIntGenerator(0, positionArray.size()-1) };
	translation.index = i;

	std::vector<double> positionParticle ( positionArray[i] );
	std::vector<double> randomVector ( randomVectorDoubleGenerator(3, -rbox, rbox) );
	positionParticle = vectorSum (positionParticle, randomVector);
	positionParticle = periodicBC (positionParticle, lengthCube);
	translation.posParticle = positionParticle;
	return translation;

}

enePos mcMove(enePos ep, std::vector<double> radiusArray, std::vector<std::vector<int>> neighborList, double rc, double lengthCube, double temp, double rbox)
{
	/*
	 *This function implements a Monte Carlo move: translation of a random particle, calculation of the energy of the new system and then acceptation or not of the move.
	 */
	indexPosParticle translation ( mcTranslation(ep.posMatrix, rbox, lengthCube) );
	std::vector<int> neighborIList = neighborList[translation.index];
	double oldEnergyParticle { energyParticle(1, translation.index, ep.posMatrix[translation.index], ep.posMatrix, neighborIList, radiusArray, rc, lengthCube) };
	double newEnergyParticle { energyParticle(0, translation.index, translation.posParticle, ep.posMatrix, neighborIList, radiusArray, rc, lengthCube) };
	bool acceptMove { metropolis (newEnergyParticle, oldEnergyParticle, temp) };

	if (acceptMove)
	{
		double newEnergy {ep.ene - oldEnergyParticle + newEnergyParticle};
		ep.ene = newEnergy;
		ep.posMatrix[translation.index] = translation.posParticle;
	}

	return ep;
}


void mcTotal(std::vector<std::vector<double>> positionArray, std::vector<double> radiusArray, double rc, double lengthCube, double temp, double rbox, std::string number)

/* This function is the core of  Monte Carlo program. It iterates over a certain amount of steps that the user can choose.
 * At each step a Monte Carlo move is made. For now the Monte Carlo move is a translation.
 * The position of each particle is written in a .xyz file at the end of each step.
  */

{

	double energy {energySystem(positionArray, radiusArray, rc, lengthCube)};

	std::vector<std::vector<int>> neighborList ( createNeighborList(positionArray, 3., lengthCube) );

	enePos ep;
	ep.ene = energy;
	ep.posMatrix = positionArray;
	std::string prename {"/home/rsimon/Documents/output/outXYZ/n600t2rho1,2/position"};
	std::string extname {".xyz"};

	for (int i = 0; i < 100000; i++)
	{
		std::cout << i <<"\n";
		ep = mcMove(ep, radiusArray, neighborList, rc, lengthCube, temp, rbox);
		if ((i % 30) == 0)
		{
			neighborList = createNeighborList(ep.posMatrix, 3., lengthCube);
		}
		if (i > 490000)
		{
			//saveInXYZ(ep.posMatrix,  radiusArray, prename + std::to_string(i) + extname );
		}

		saveEnergyTXT(ep.ene, "/Users/romainsimon/eclipse-workspace/swapMC/output/outE/quicktest.txt");

	}
}
