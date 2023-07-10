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
    std::string energyFilePath{m_folderPath + "/outE.txt"};
    std::string pressureFilePath {m_folderPath + "/outP.txt"};
    // std::vector<double> radiusArray (divideVectorByScalar(m_typeArray, 2));

	//saveInXYZ(m_positionArray, radiusArray, m_moleculeType, m_lengthCube,preName + std::to_string(0) + extname );
	saveDoubleTXT(m_energy / m_nParticles, m_folderPath + "/outE.txt");
	// saveDisplacement(m_totalDisplacementMatrix, preNameDisp + std::to_string(0) + extnameDisp);

	const std::vector<int> saveTimeStepArray ( createSaveTime(m_timeSteps, m_saveUpdate, 1.1));

	int save_index { 0 };


	for (int i = 0; i < m_timeSteps; i++) //Iteration over m_timeSteps
	{

		for (int j = 0; j < m_nParticles; j++) // N Monte Carlo moves are tried in one time step.
		{
			mcMove();
		}

		//m_interDisplacementMatrix = matrixSum( m_interDisplacementMatrix, m_stepDisplacementMatrix );
		//m_totalDisplacementMatrix = matrixSum( m_totalDisplacementMatrix, m_stepDisplacementMatrix );
		//std::fill( m_stepDisplacementMatrix.begin(), m_stepDisplacementMatrix.end(), std::vector<double>(3, 0));

		if (m_neighMethod == "verlet")
		{
            if (m_simulationMol == "polymer")
            {
                m_systemNeighbors.checkInterDisplacement(m_systemParticles,m_bondPotentials);
            }
            else
            {
                m_systemNeighbors.checkInterDisplacement(m_systemParticles);
            }
		}

		// Next, the results of the simulations are saved.

		if (saveTimeStepArray[save_index] == i)
		{
			//radiusArray = divideVectorByScalar(m_typeArray, 2);
            std::string nameXYZ {preName};
            nameXYZ.append(std::to_string(i + 1)).append(extname);
            std::string nameDisp {preNameDisp};
            nameDisp.append(std::to_string(i + 1)).append(extnameDisp);
            //saveInXYZ (m_positionArray, radiusArray, m_moleculeType, m_lengthCube, nameXYZ);
			//saveDisplacement (m_totalDisplacementMatrix, nameDisp);
			++save_index;
		}

		if (i % m_saveRate == 0)
		{
			saveDoubleTXT(m_energy / m_nParticles, energyFilePath); //Energy is saved at each time step.
		}
		if (m_calculatePressure) // Saves pressure if m_calculatePressure is True
		{
			saveDoubleTXT( m_pressure, pressureFilePath);
		}

	}

    //radiusArray = divideVectorByScalar(m_typeArray, 2);
	m_acceptanceRate /= m_timeSteps;
    //saveInXYZ(m_positionArray, radiusArray, m_moleculeType, m_lengthCube, preName + std::to_string(m_timeSteps) + extname );
    //saveDisplacement(m_totalDisplacementMatrix, preNameDisp + std::to_string(m_timeSteps) + extnameDisp);
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
												   m_positionArray, neighborIList, m_typeArray, bondsI,
												   m_squareRc, m_lengthCube, m_squareR0, m_feneK);
		newEnergyParticle = energyParticlePolymer (indexTranslation, tentativePositionArray[indexTranslation + j],
												   tentativePositionArray, neighborIList, m_typeArray, bondsI,
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
    std::vector<double> positionParticleTranslation = m_systemParticles.getPositionI(indexTranslation);
    positionParticleTranslation = vectorSum(positionParticleTranslation, randomVector);
    positionParticleTranslation = periodicBC(positionParticleTranslation, m_systemParticles.getLengthCube());

    std::vector<int> neighborIList { m_systemNeighbors.getNeighborIList(indexTranslation) };
    double oldEnergyParticle;
    double newEnergyParticle;

    if (m_simulationMol == "polymer")
    {

        oldEnergyParticle = energyParticlePolymer(indexTranslation,
                                                  m_systemParticles.getPositionI(indexTranslation),
                                                  m_systemParticles, neighborIList,
                                                  m_pairPotentials, m_bondPotentials);
        newEnergyParticle = energyParticlePolymer(indexTranslation, positionParticleTranslation,
                                                  m_systemParticles, neighborIList,
                                                  m_pairPotentials, m_bondPotentials);

    }
    else
    {
        oldEnergyParticle = energyParticle(indexTranslation,
                                           m_systemParticles.getPositionI(indexTranslation),
                                           m_systemParticles, neighborIList,
                                           m_pairPotentials);

        newEnergyParticle = energyParticle(indexTranslation, positionParticleTranslation,
                                           m_systemParticles, neighborIList,m_pairPotentials);
    }

    const double diff_energy {newEnergyParticle - oldEnergyParticle};

    // Metropolis criterion
    const bool acceptMove{ metropolis(diff_energy) };

    // If the move is accepted, then the energy, the position array and the displacement array can be updated.
    // If the m_calculatePressure is set to True, then the pressure is calculated.
    if (acceptMove)
    {
        generalUpdate(diff_energy);
        /***
        if (m_calculatePressure)
        {
            const double newPressureParticle {pressureParticle(m_temp, indexTranslation, positionParticleTranslation, m_positionArray, neighborIList, m_typeArray, m_squareRc, m_lengthCube, m_halfLengthCube)};
            const double oldPressureParticle {pressureParticle(m_temp, indexTranslation, m_positionArray[indexTranslation], m_positionArray, neighborIList, m_typeArray, m_squareRc, m_lengthCube, m_halfLengthCube)};
            m_pressure += newPressureParticle - oldPressureParticle;
        }
        ***/
        m_systemNeighbors.updateInterDisplacement(indexTranslation, randomVector);
        m_systemParticles.updatePositionI(indexTranslation, positionParticleTranslation);
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

    double energyParticle1;
    double energyParticle2;
    double energyParticleSwap1;
    double energyParticleSwap2;
    const std::vector<int> neighborIList1 { m_systemNeighbors.getNeighborIList(indexSwap1) };
    const std::vector<int> neighborIList2 { m_systemNeighbors.getNeighborIList(indexSwap2) };


    if (m_simulationMol == "polymer")
    {

        energyParticle1 = energyParticlePolymer(indexSwap1,m_systemParticles.getPositionI(indexSwap1),
                                                m_systemParticles, neighborIList1,
                                                m_pairPotentials, m_bondPotentials, indexSwap2);

        energyParticle2 = energyParticlePolymer(indexSwap2, m_systemParticles.getPositionI(indexSwap2),
                                                m_systemParticles, neighborIList2,
                                                m_pairPotentials, m_bondPotentials, indexSwap1);

        m_systemParticles.swapParticleTypesIJ(indexSwap1, indexSwap2);

        energyParticleSwap1 = energyParticlePolymer(indexSwap1,m_systemParticles.getPositionI(indexSwap1),
                                                    m_systemParticles, neighborIList1,
                                                    m_pairPotentials, m_bondPotentials, indexSwap2);

        energyParticleSwap2 = energyParticlePolymer(indexSwap2, m_systemParticles.getPositionI(indexSwap2),
                                                    m_systemParticles, neighborIList2,
                                                    m_pairPotentials, m_bondPotentials, indexSwap1);

    }
    else
    {
        energyParticle1 = energyParticle(indexSwap1, m_systemParticles.getPositionI(indexSwap1),
                                         m_systemParticles, neighborIList1,
                                         m_pairPotentials, indexSwap2);

        energyParticle2 = energyParticle(indexSwap2, m_systemParticles.getPositionI(indexSwap2),
                                         m_systemParticles, neighborIList2,
                                         m_pairPotentials, indexSwap1);

        m_systemParticles.swapParticleTypesIJ(indexSwap1, indexSwap2);

        energyParticleSwap1 = energyParticle(indexSwap1, m_systemParticles.getPositionI(indexSwap1),
                                             m_systemParticles, neighborIList1,
                                             m_pairPotentials, indexSwap2);

        energyParticleSwap2 = energyParticle(indexSwap2, m_systemParticles.getPositionI(indexSwap2),
                                             m_systemParticles, neighborIList2,
                                             m_pairPotentials, indexSwap1);
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
        /***
        if (m_calculatePressure)
        {
            const double pressureParticleSwap1 {pressureParticle(m_temp, indexSwap1,
                                                           m_positionArray[indexSwap1],
                                                           m_positionArray, neighborIList1,
                                                           m_typeArray, m_squareRc,
                                                           m_lengthCube, m_halfLengthCube)};

            const double pressureParticleSwap2 {pressureParticle(m_temp, indexSwap2,
                                                           m_positionArray[indexSwap2],
                                                           m_positionArray, neighborIList2,
                                                           m_typeArray, m_squareRc,
                                                           m_lengthCube, m_halfLengthCube)};

            m_typeArray[indexSwap1] = type1;
            m_typeArray[indexSwap2] = type2;

            const double pressureParticle1 {pressureParticle(m_temp, indexSwap1,
                                                           m_positionArray[indexSwap1],
                                                           m_positionArray, neighborIList1,
                                                           m_typeArray, m_squareRc,
                                                           m_lengthCube, m_halfLengthCube)};


            const double pressureParticle2 {pressureParticle(m_temp, indexSwap2,
                                                           m_positionArray[indexSwap2],
                                                           m_positionArray, neighborIList2,
                                                           m_typeArray, m_squareRc,
                                                           m_lengthCube, m_halfLengthCube)};

            m_typeArray[indexSwap1] = type2;
            m_typeArray[indexSwap2] = type1;

            m_pressure += pressureParticleSwap1 + pressureParticleSwap2 - pressureParticle1 - pressureParticle2;
        }
        ***/
    }
    else
    {
        m_systemParticles.swapParticleTypesIJ(indexSwap1, indexSwap2);
    }
}


void MonteCarlo::generalUpdate(double diff_energy)
{

    const double newEnergy {m_energy + diff_energy};
    m_energy = newEnergy;
    m_acceptanceRate += 1. / m_nParticles; // increment of the acceptance rate.

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