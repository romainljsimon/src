/*
 * energy.cpp
 *
 *  Created on: 3 oct. 2022
 *      Author: romainsimon
 */
#include <iostream>
#include <vector>
#include <cmath>
#include "util.h"
#include <limits>

/*******************************************************************************
 * This function calculates the Lennard-Jones potential energy between two
 * particles seperated by a distance whose square is equal to squareDistance.
 * We consider a cut off whose square is equal to squareRc times the mean diameter
 * of the two particles considered. The potential is equal to 0 for distances
 * greater that the cut off.
 *
 * @param squareDistance Square of the distance seperating the two considered
 *                       particles
 *        sigmaA Particle A's diameter
 *        sigmaB Particle B's diameter
 *        squareRc Square of the cut off radius
 *        shift Shifts the Lennard-Jones potential by a constant value
 *
 * @return energy
 ******************************************************************************/
double ljPotential(const double& squareDistance, const double& sigmaA, const double& sigmaB, const double& squareRc, const double& shift)
{

	double squareSigma { std::pow (((sigmaA + sigmaB) / 2. ), 2.) };
    // We consider a cut-off radius which is the threshold maximum distance of interaction between two particles
	if (squareDistance > squareRc * squareSigma)
	{
		return 0.;
	}

	else
	{
		double rap_square {std::pow ((squareSigma / squareDistance), 3.)};
		return 4. * rap_square * ( rap_square - 1.) + 4. * shift;
	}
}

/*******************************************************************************
 * This function calculates the FENE (finitely extensible nonlinear elastic) 
 * potential energy between two nearest-neighbor monomers seperated by a 
 * distance whose square is equal to squareDistance. SquareR0 and feneK are two
 * parameters of the potential. SquareRo is the maximum square distance that can
 * seperate two neighboring atoms and feneK represents the stiffness of the 
 * elastic. They are chosen such as bond crossing doesn't occur.
 *
 * @param squareDistance Square of the distance seperating the two considered
 *                       particles
 *        sigmaA particle A's diameter.
 *        sigmaB particle B's diameter.
 *        squareR0 maximum square distance.
 *        feneK stiffness of the elastic.
 *
 * @return energy
 ******************************************************************************/
double fenePotential(const double& squareDistance, const double& sigmaA, const double& sigmaB, double squareR0, double feneK)
{
	double squareSigma { std::pow (((sigmaA + sigmaB) / 2. ), 2.) };
	squareR0 = squareR0 * squareSigma;
    // We consider a cut-off radius which is the threshold maximum distance of interaction between two particles
    if (squareDistance >= squareR0)
	{
		return std::numeric_limits<double>::infinity();
	}

	else
	{
		feneK = feneK / squareSigma;
		return -0.5 * feneK * squareR0 * std::log(1. - squareDistance / squareR0);
	}
}

/*******************************************************************************
 * This function calculates the distance between two particles considering the 
 * periodic boundary conditions.
 *
 * @param positionA particle A's position.
 *        particleB particle B's position.
 *        lengthCube simulation box length.
 *
 * @return distance between two particles.
 ******************************************************************************/
double squareDistancePair(const std::vector<double>& positionA,  const std::vector<double>& positionB, const double& lengthCube)
{
	double squareDistance { 0. };
	double halfLengthCube { lengthCube / 2. };

	for (int i = 0; i < 3; i++)
	{
		double diff { positionA[i] - positionB[i] };

		if (diff > halfLengthCube)
			diff -= lengthCube;

		else if (diff < - halfLengthCube)
			diff += lengthCube;

		squareDistance += std::pow(diff, 2);
	}

	return squareDistance;
}


/*******************************************************************************
 * This function calculates the potential energy of one particle considering
 * that particles are Lennard-Jones particles.
 *
 * @param indexParticle Considered particle's index in the positionArray and 
 *                      diameterArray arrays.
 *        positionParticle Considered particle's position.
 * 	      positionArray array of particle positions.
 * 	      neighborIList Considered particle's neighbor list.
 *        diameterArray array of particle diameters.
 *        squareRc square of the cut off radius.
 *        lengthCube simulation box length.
 * @return Particle's total energy.
 ******************************************************************************/
double energyParticle(const int& indexParticle, const std::vector<double>& positionParticle,
		      const std::vector<std::vector<double>>& positionArray, const std::vector<int>& neighborIList,
              const std::vector<double>& diameterArray, const double& squareRc, const double& lengthCube)
{
	double energy { 0. };
	double particleDiameter = diameterArray[indexParticle];
	int neighborIListSize { static_cast<int>(neighborIList.size()) };

	for (int i = 0; i < neighborIListSize; i++)
	{
		int realIndex = neighborIList[i];
		double squareDistance { squareDistancePair (positionParticle, positionArray[realIndex], lengthCube)};

		if (realIndex == indexParticle)
		{
			continue;
		}

		energy += ljPotential(squareDistance, particleDiameter, diameterArray[realIndex], squareRc, 0.);
	}
	return energy;
}


/*******************************************************************************
 * This function calculates the total potential energy of the system considering
 * that particles are Lennard-Jones particles.
 *
 * @param positionArray array of particle positions.
 *        diameterArray array of particle diameter.
 *        squareRc square of the cut off radius.
 *        lengthCube simulation box length.
 *
 * @return System's total energy.
 ******************************************************************************/
double energySystem(const std::vector<std::vector<double>>& positionArray, const std::vector<double>& diameterArray,
		            const std::vector<std::vector<int>>& neighborList, const double& squareRc,
					const double& lengthCube)
{
	double energy { 0. };
	int positionArraySize {static_cast<int>( positionArray.size() )};

	for (int i = 0; i < positionArraySize - 1; i++) //Outer loop for rows
    {
    	std::vector<int> neighborIList ( neighborList[i] );

    	energy += energyParticle (i, positionArray[i], positionArray, neighborIList,
				  diameterArray, squareRc, lengthCube) / 2;

    	/***
		for (int j = i + 1; j < positionArraySize; j++) //inner loop for columns
        {
			double squareDistance { squareDistancePair (positionArray[i], positionArray[j], lengthCube)};
        	energy += ljPotential(squareDistance, 2. * diameterArray[j], 2. * diameterArray[i], squareRc, 0.);
        }
        ***/
    }
	return energy;
}

/*******************************************************************************
 * This function calculates the potential energy of one monomer considering
 * that it is part of a polymer chain. Every monomer interact via a 
 * Lennard-Jones potential. Nearest-neighbor monomers also interact via a FENE 
 * potential.
 *
 * @param indexParticle Considered monomer's index in the positionArray and 
 *                      diameterArray arrays.
 *        positionArray Array of monomer positions.
 *        neighborIList Considered monomer's neighbor list.
 *        diameterArray Array of particle radiuses.
 *        bondsI vector that indicates which monomers are linked with the 
 *               considered monomer.
 *        squareRc Square of the cut off radius.
 *        lengthCube Simulation box length.
 *        squareR0 Maximum square distance.
 *        feneK Elastic's stiffness.
 *
 * @return Polymeric particle's total energy.
 ******************************************************************************/
double energyParticlePolymer (const int& indexParticle, const std::vector<double>& positionParticle,
							  const std::vector<std::vector<double>>& positionArray, const std::vector<int>& neighborIList,
							  const std::vector<double>& diameterArray, const std::vector<int>& bondsI,
							  const double& squareRc, const double& lengthCube, const double& squareR0, const double& feneK)
{
	double energy { 0. };
	double particleDiameter = diameterArray[indexParticle];
	int neighborIListSize { static_cast<int>(neighborIList.size()) };
	int bondsISize { static_cast<int>(bondsI.size()) };

	for (int i = 0; i < neighborIListSize; i++)
	{
		int realIndex = neighborIList[i];
		double squareDistance { squareDistancePair (positionParticle, positionArray[realIndex], lengthCube)};

		if (realIndex == indexParticle)
		{
			continue;
		}

		energy += ljPotential(squareDistance, particleDiameter, diameterArray[realIndex], squareRc, 0.25); //127. / 4096);
	}

	for (int i = 0; i < bondsISize; i++)
	{
		int realIndex = bondsI[i];

		if ((realIndex == indexParticle) || (realIndex == -1))
		{
			continue;
		}

		/*** Weak Fene case
		if (((( indexParticle % 3 ) == 0) && (realIndex == indexParticle + 2)) || (( indexParticle == realIndex + 2 ) && ( (realIndex % 3 ) == 0)))
		{

			double squareDistance { squareDistancePair (positionParticle, positionArray[realIndex], lengthCube)};
			energy += fenePotential(squareDistance, 2. * particleRadius, 2. * diameterArray[realIndex], 1.56, 0.3);
		}
		else
		{
			double squareDistance { squareDistancePair (positionParticle, positionArray[realIndex], lengthCube)};
			energy += fenePotential(squareDistance, 2. * particleRadius, 2. * diameterArray[realIndex], squareR0, feneK);
		}
		***/

		/*** Normal case


		double squareDistance { squareDistancePair (positionParticle, positionArray[realIndex], lengthCube)};
		energy += fenePotential(squareDistance, 2. * particleRadius, 2. * diameterArray[realIndex], squareR0, feneK);
		***/

		/*** Shifted WCA model

		double squareDistance { squareDistancePair (positionParticle, positionArray[realIndex], lengthCube)};

		if (((( indexParticle % 3 ) == 0) && (realIndex == indexParticle + 2)) || (( indexParticle == realIndex + 2 ) && ( (realIndex % 3 ) == 0)))
		{

			energy -= ljPotential(squareDistance, 2. * particleRadius, 2. * diameterArray[realIndex], squareRc, 0.25);
			energy += ljPotential(squareDistance - 0.25, 2. * particleRadius, 2. * diameterArray[realIndex], squareRc, 0.25);
		}

		energy += fenePotential(squareDistance, 2. * particleRadius, 2. * diameterArray[realIndex], squareR0, feneK);
		***/

		/*** Good shifted model

		double squareDistance { squareDistancePair (positionParticle, positionArray[realIndex], lengthCube)};

		if (((( indexParticle % 3 ) == 0) && (realIndex == indexParticle + 2)) || (( indexParticle == realIndex + 2 ) && ( (realIndex % 3 ) == 0)))
		{

			energy -= ljPotential(squareDistance, 2. * particleRadius, 2. * diameterArray[realIndex], squareRc, 0.25);
			energy += ljPotential(squareDistance - 0.4, 2. * particleRadius, 2. * diameterArray[realIndex], squareRc, 0.25);
			energy += fenePotential(squareDistance, 2. * particleRadius, 2. * diameterArray[realIndex], 2.89, feneK);
		}
		else
		{
			energy += fenePotential(squareDistance, 2. * particleRadius, 2. * diameterArray[realIndex], squareR0, feneK);
		}
		***/

		// Large End case
		double squareDistance { squareDistancePair (positionParticle, positionArray[realIndex], lengthCube)};


		if (((( indexParticle % 3 ) == 0) && (realIndex == indexParticle + 2)) || (( indexParticle == realIndex + 2 ) && ( (realIndex % 3 ) == 0)))
		{
			double sigma_factor = 1.35;
			double newSquareR0 = 1;
			energy -= ljPotential(squareDistance, particleDiameter, diameterArray[realIndex], squareRc, 0.25);
			energy += ljPotential(squareDistance, particleDiameter * sigma_factor, diameterArray[realIndex] * sigma_factor, squareRc, 0.25);
			energy += fenePotential(squareDistance, particleDiameter * sigma_factor, diameterArray[realIndex] * sigma_factor, squareR0, feneK);
		}
		else
		{
			energy += fenePotential(squareDistance, particleDiameter , diameterArray[realIndex], squareR0, feneK);
		}

	}

	return energy;
}


/*******************************************************************************
 * This function calculates the total potential energy of the system considering
 * that the system is polymeric. Every monomer interact via a Lennard-Jones
 * potential. Nearest-neighbor monomers also interact via a FENE potential.
 *
 * @param positionArray array of particle positions.
 *        diameterArray array of particle diameters.
 *        bondsMatrix matrix that indicates which monomers are linked together.
 *        squareRc square of the cut off radius.
 *        lengthCube simulation box length.
 *        squareR0 maximum square distance.
 *        feneK stiffness of the elastic
 *
 * @return Polymeric system's total energy.
 ***********************************************************position65000.xyz*******************/
double energySystemPolymer(const std::vector<std::vector<double>>& positionArray, const std::vector<double>& diameterArray,
						   const std::vector<std::vector<int>>& bondsMatrix, const std::vector<std::vector<int>>& neighborList,
						   const double& squareRc, const double& lengthCube, const double& squareR0, const double& feneK)
{
	double energy { 0. };
	int positionArraySize {static_cast<int>(positionArray.size())};

	for (int i = 0; i < positionArraySize; i++) //Outer loop for rows
    {

    	std::vector<int> bondsI ( bondsMatrix[i] );
    	int bondsISize { static_cast<int>(bondsI.size()) };

    	std::vector<int> neighborIList ( neighborList[i] );

    	energy += energyParticlePolymer (i, positionArray[i], positionArray, neighborIList,
				  diameterArray, bondsI, squareRc, lengthCube, squareR0, feneK) / 2.;

    	/***
		for (int j = i + 1; j < positionArraySize; j++) //inner loop for columns
        {
			double squareDistance { squareDistancePair (positionArray[i], positionArray[j], lengthCube)};

			energy += ljPotential(squareDistance, 2 * diameterArray[j], 2 * diameterArray[i], squareRc, 0.25); //127. / 4096);

        }
		for (int j = 0; j < bondsISize; j++)
		{
			int realIndex = bondsI[j];

			if ((realIndex == -1) || (realIndex < i))
			{
				continue;
			}

			double squareDistance { squareDistancePair (positionArray[i], positionArray[realIndex], lengthCube)};
			energy += fenePotential(squareDistance, 2. * diameterArray[i], 2. * diameterArray[realIndex], squareR0, feneK);

        }
        ***/
    }
	return energy;
}

double rforce(double squareDistance, double sigmaA, double sigmaB, double squareRc)
{
	double squareSigma { std::pow (((sigmaA + sigmaB) / 2. ), 2.) };
	if (squareDistance > squareRc * squareSigma) // We consider a cut-off radius which is the threshold maximum distance of interaction between two particles
	{
		return 0.;
	}
	else
	{
		return 24. * (std::pow (( squareSigma / squareDistance), 3. ) - 2. * std::pow (( squareSigma / squareDistance ), 6.));
	}
}


double pressureParticle(double temp, int indexParticle, std::vector<double> positionParticle,
		std::vector<std::vector<double>> positionArray, std::vector<int> neighborIList, std::vector<double> diameterArray,
		double squareRc, double lengthCube)
{
	double pressure { 0. };
	double particleDiameter = diameterArray[indexParticle];
	int neighborIListSize { static_cast<int>(neighborIList.size()) };
	int nParticles {static_cast<int>(positionArray.size())};

	for (int i = 0; i < neighborIListSize; i++)
	{
		int realIndex = neighborIList[i];
		double squareDistance { squareDistancePair (positionParticle, positionArray[realIndex], lengthCube) };

		if (realIndex == indexParticle)
		{
			continue;
		}

		pressure += rforce(squareDistance, particleDiameter, diameterArray[realIndex], squareRc);
	}
	return  - pressure / (3. * temp *  nParticles);
}


double pressureSystem(double temp, std::vector<std::vector<double>> positionArray, std::vector<double> diameterArray, double squareRc, double lengthCube)
{
	double sum { 0 };
	int positionArraySize {static_cast<int>(positionArray.size())};
	for (int i = 0; i < positionArraySize - 1; i++)
	{
		for (int j = i + 1; j < positionArraySize; j++) //inner loop for columns
        {
			double squareDistance { squareDistancePair (positionArray[i], positionArray[j], lengthCube)};
			sum += rforce(squareDistance, 2. * diameterArray[i], 2. * diameterArray[j], squareRc);
        }

	}
	return 1. - sum / (3. *  positionArraySize * temp);
}

