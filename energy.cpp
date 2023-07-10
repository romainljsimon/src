/*
 * energy.cpp
 *
 *  Created on: 3 oct. 2022
 *      Author: Romain Simon
 */
#include <vector>
#include "POTENTIALS/PairPotentials.h"
#include "POTENTIALS/BondPotentials.h"
#include "PARTICLES/Particles.h"
#include "NEIGHBORS/Neighbors.h"

/*******************************************************************************
 * This function calculates the Lennard-Jones potential energy between two
 * particles seperated by a distance whose square is equal to squareDistance.
 * We consider a cut off whose square is equal to squareRc times the mean diameter
 * of the two particles considered. The potential is equal to 0 for distances
 * greater that the cut off.
 *
 * @param squareDistance Square of the distance separating the two considered
 *                       particles
 *        sigmaA Particle A's diameter
 *        sigmaB Particle B's diameter
 *        squareRc Square of the cut off radius
 *        shift Shifts the Lennard-Jones potential by a constant value
 *
 * @return energy
 ******************************************************************************/


/*******************************************************************************
 * This function calculates the FENE (finitely extensible nonlinear elastic) 
 * potential energy between two nearest-neighbor monomers seperated by a 
 * distance whose square is equal to squareDistance. SquareR0 and feneK are two
 * parameters of the potential. SquareRo is the maximum square distance that can
 * separate two neighboring atoms and feneK represents the stiffness of the
 * elastic. They are chosen such as bond crossing doesn't occur.
 *
 * @param squareDistance Square of the distance separating the two considered
 *                       particles
 *        sigmaA particle A's diameter.
 *        sigmaB particle B's diameter.
 *        squareR0 maximum square distance.
 *        feneK stiffness of the elastic.
 *
 * @return energy
 ******************************************************************************/

double energyParticle(const int& indexParticle, const std::vector<double>& positionParticle,
                      const Particles& systemParticles, const std::vector<int>& neighborIList,
                      const PairPotentials& pairPotentials, const int& indexSkip)
{
    return pairPotentials.energyPairParticle(indexParticle, positionParticle, systemParticles, neighborIList, indexSkip);
}

double energyParticle(const int& indexParticle, const std::vector<double>& positionParticle,
                      const Particles& systemParticles, const std::vector<int>& neighborIList,
                      const PairPotentials& pairPotentials)
{
    return pairPotentials.energyPairParticle(indexParticle, positionParticle, systemParticles, neighborIList, -1);
}


/*******************************************************************************
 * This function calculates the total potential energy of the system considering
 * that particles are Lennard-Jones particles.
 *
 * @param positionArray array of particle positions.
 *        typeArray array of particle diameter.
 *        squareRc square of the cut off radius.
 *        lengthCube simulation box length.
 *
 * @return System's total energy.
 ******************************************************************************/
double energySystem(const Particles& systemParticles, const Neighbors& systemNeighbors,
                           const PairPotentials& pairPotentials)
{
    double energy { 0. };
    const int positionArraySize { systemParticles.getNParticles() };

    for (int indexParticle = 0; indexParticle < positionArraySize; indexParticle++) //Outer loop for rows
    {
        const std::vector<int>& neighborIList (systemNeighbors.getNeighborIList(indexParticle) );
        energy += energyParticle (indexParticle, systemParticles.getPositionI(indexParticle),
                                         systemParticles, neighborIList, pairPotentials) / 2.;
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
 *                      typeArray arrays.
 *        positionArray Array of monomer positions.
 *        neighborIList Considered monomer's neighbor list.
 *        typeArray Array of particle radii.
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
							  const Particles& systemParticles, const std::vector<int>& neighborIList,
							  const PairPotentials& pairPotentials, const BondPotentials& bondPotentials,
                              const int& indexSkip)
{
	double energy { 0. };
    energy += pairPotentials.energyPairParticle(indexParticle, positionParticle, systemParticles,
                                                neighborIList, indexSkip);

    energy += bondPotentials.feneBondEnergyI(indexParticle, positionParticle, systemParticles, indexSkip);
    return energy;
}

double energyParticlePolymer (const int& indexParticle, const std::vector<double>& positionParticle,
                              const Particles& systemParticles, const std::vector<int>& neighborIList,
                              const PairPotentials& pairPotentials, const BondPotentials& bondPotentials)
{
    double energy { 0. };
    energy += pairPotentials.energyPairParticle(indexParticle, positionParticle, systemParticles,
                                                neighborIList, -1);
    energy += bondPotentials.feneBondEnergyI(indexParticle, positionParticle, systemParticles, -1);
    return energy;
}


/*******************************************************************************
 * This function calculates the total potential energy of the system considering
 * that the system is polymeric. Every monomer interact via a Lennard-Jones
 * potential. Nearest-neighbor monomers also interact via a FENE potential.
 *
 * @param positionArray array of particle positions.
 *        typeArray array of particle diameters.
 *        bondsMatrix matrix that indicates which monomers are linked together.
 *        squareRc square of the cut off radius.
 *        lengthCube simulation box length.
 *        squareR0 maximum square distance.
 *        feneK stiffness of the elastic
 *
 * @return Polymeric system's total energy.
 ***********************************************************position65000.xyz*******************/
double energySystemPolymer(const Particles& systemParticles, const Neighbors& systemNeighbors,
                           const PairPotentials& pairPotentials, const BondPotentials& bondPotentials)
{
	double energy { 0. };
	const int positionArraySize { systemParticles.getNParticles() };
	for (int indexParticle = 0; indexParticle < positionArraySize; indexParticle++) //Outer loop for rows
    {
    	const std::vector<int>& neighborIList (systemNeighbors.getNeighborIList(indexParticle) );
    	energy += energyParticlePolymer (indexParticle, systemParticles.getPositionI(indexParticle),
                                         systemParticles, neighborIList, pairPotentials, bondPotentials) / 2.;
        //std::cout << energy << "\n";
    }
	return energy;
}
