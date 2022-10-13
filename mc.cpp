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

struct enePos
{
	double ene;
	std::vector<std::vector<double>> posMatrix;
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

double periodicBC(double position, double lengthCube)
/*
 *This function is an implementation of the periodic Boundary conditions.
 *If a particle gets out of the simulation box from one of the sides, it gets back in the box from the opposite side.
 */
{

	if (position < 0)
		position += lengthCube;

	else if (position > lengthCube)
		position -= lengthCube;

	return position;
}

std::vector<std::vector <double>> mcTranslation(std::vector<std::vector<double>>& positionArray, double rbox, double lengthCube)
{
	/*
	 *This function chooses a particle randomly and translates it randomly inside of box of size rbox.
	 */
	int i {randomIntGenerator(0, positionArray.size()-1)};

	std::vector<std::vector <double>> newPositionArray (positionArray);

	double newPositionX = newPositionArray[i][0] + randomDoubleGenerator(-rbox, rbox);
	double newPositionY = newPositionArray[i][1] + randomDoubleGenerator(-rbox, rbox);
	double newPositionZ = newPositionArray[i][2] + randomDoubleGenerator(-rbox, rbox);

	newPositionArray[i][0] = periodicBC(newPositionX, lengthCube);
	newPositionArray[i][1] = periodicBC(newPositionY, lengthCube);
	newPositionArray[i][2] = periodicBC(newPositionZ, lengthCube);

	return newPositionArray;

}

enePos mcMove(enePos ep, std::vector<double> radiusArray, double rc, double lengthCube, double temp, double rbox)
{
	/*
	 *This function implements a Monte Carlo move: translation of a random particle, calculation of the energy of the new system and then acceptation or not of the move.
	 */
	std::vector<std::vector <double>> newPositionArray ( mcTranslation(ep.posMatrix, rbox, lengthCube));
	double newEnergy {energySystem(newPositionArray, radiusArray, rc, lengthCube)};
	bool acceptMove {metropolis (newEnergy, ep.ene, temp)};

	if (acceptMove)
	{
		ep.ene = newEnergy;
		ep.posMatrix = newPositionArray;
	}

	return ep;
}


void mcTotal(std::vector<std::vector<double>> positionArray, std::vector<double> radiusArray, double rc, double lengthCube, double temp, double rbox)

/* This function is the core of  Monte Carlo program. It iterates over a certain amount of steps that the user can choose.
 * At each step a Monte Carlo move is made. For now the Monte Carlo move is a translation.
 * The position of each particle is written in a .xyz file at the end of each step.
  */

{
	double energy {energySystem(positionArray, radiusArray, rc, lengthCube)};
	enePos ep;
	ep.ene = energy;
	ep.posMatrix = positionArray;
	std::string prename {"output/outXYZ/n300rho1t1v2/position"};
	std::string extname {".xyz"};

	for (int i = 0; i < 10000; i++)
	{
		ep = mcMove(ep, radiusArray, rc, lengthCube, temp, rbox);
		saveInXYZ(ep.posMatrix,  radiusArray, prename + std::to_string(i) + extname );
		saveEnergyTXT(ep.ene, "output/outE/n300rho1t1v2.txt");
	}
}
