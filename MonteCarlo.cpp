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

	const std::string preName ("./outXYZ/position");
	const std::string extname {".xyz"};
	const std::string preNameDisp ("./disp/displacement");
	const std::string extnameDisp {".txt"};
    constexpr std::string_view energyFilePath{"./outE.txt"};
    //std::string pressureFilePath {"./outP.txt"};
    // std::vector<double> radiusArray (divideVectorByScalar(m_typeArray, 2));

	m_systemMolecules.saveInXYZ(preName + std::to_string(0) + extname );
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

        m_systemNeighbors.checkInterDisplacement(m_systemMolecules);


		// Next, the results of the simulations are saved.

		if (saveTimeStepArray[save_index] == i)
		{
			//radiusArray = divideVectorByScalar(m_typeArray, 2);
            std::string nameXYZ {preName};
            nameXYZ.append(std::to_string(i + 1)).append(extname);
            //std::string nameDisp {preNameDisp};
            //nameDisp.append(std::to_string(i + 1)).append(extnameDisp);
            m_systemMolecules.saveInXYZ(nameXYZ );
			//saveDisplacement (m_totalDisplacementMatrix, nameDisp);
			++save_index;
		}

		if (i % m_saveRate == 0)
		{
			saveDoubleTXT(m_energy / m_nParticles, energyFilePath); //Energy is saved at each time step.
		}
        /***
		if (m_calculatePressure) // Saves pressure if m_calculatePressure is True
		{
			saveDoubleTXT( m_pressure, pressureFilePath);
		}
        ***/

	}

    //radiusArray = divideVectorByScalar(m_typeArray, 2);
	m_acceptanceRate /= m_timeSteps;
    m_systemMolecules.saveInXYZ( preName + std::to_string(m_timeSteps) + extname );
    //saveDisplacement(m_totalDisplacementMatrix, preNameDisp + std::to_string(m_timeSteps) + extnameDisp);
    // saveDoubleTXT( m_errors, m_folderPath + "/errors.txt");

    constexpr std::string_view translationString{"Translation MC move acceptance rate: "};
    if ( m_swap )
    {
        m_acceptanceRateSwap /= m_timeSteps;
        m_acceptanceRateSwap /= m_pSwap;

        const double translationRate { (m_acceptanceRate - m_pSwap * m_acceptanceRateSwap) / (1 - m_pSwap)};

        constexpr std::string_view swapString { "Swap MC move acceptance rate: " };
        constexpr std::string_view totalString { "Total MC move acceptance rate: " };

        std::cout << translationString << translationRate << "\n";
        std::cout << swapString << m_acceptanceRateSwap << "\n";
        std::cout << totalString  << m_acceptanceRate << "\n";
    }
    else
    {
        std::cout << translationString  << m_acceptanceRate << "\n";
    }
    double updateRate { static_cast<double>(m_systemNeighbors.getUpdateRate()) / m_timeSteps};
    constexpr std::string_view neighborString { "Neighbor list update rate: "};
    constexpr std::string_view neighborErrorString  { "Number of neighbor list errors: "};
	std::cout << neighborString  << updateRate << "\n";
	std::cout << neighborErrorString <<  m_systemNeighbors.getErrors() << "\n";

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
            //mcMoleculeTranslation();
            mcTranslation ();

        }
    }
    else
    {
    }
}

/*******************************************************************************
 * This function implements a Monte Carlo move: translation of a the whole
 * polymer, calculation of the energy of the new system and then acceptation or not of
 * the move according to the Metropolis criterion. It is called ??? times in one
 * time step.
 ******************************************************************************/

void MonteCarlo::mcMoleculeTranslation()
{
	const int typeMolecule { randomIntGenerator(0, (m_nParticles - 1) / 3) }; // randomly chosen Molecule
    const int indexTranslation { 3 * typeMolecule};
	std::vector<double> randomVector ( randomVectorDoubleGenerator(3, -m_rBox, m_rBox) );
	double oldEnergyMolecule {0};
	double newEnergyMolecule {0};
    std::vector<double> positionArrayTranslation;

	for (int j = 0; j < 3; j++)
	{

        const int newIndexTranslation {indexTranslation + j };
        const std::vector<double> positionParticleTranslation(vectorTranslation(newIndexTranslation,
                                                                                randomVector));
        positionArrayTranslation.insert( positionArrayTranslation.end(),
                                         positionParticleTranslation.begin(),
                                         positionParticleTranslation.end() );
        const std::vector<int>& neighborIList { m_systemNeighbors.getNeighborIList(newIndexTranslation) };

		oldEnergyMolecule += m_systemMolecules.energyPairParticleExtraMolecule( newIndexTranslation, neighborIList, typeMolecule);
		newEnergyMolecule += m_systemMolecules.energyPairParticleExtraMolecule( newIndexTranslation,
                                                                           positionParticleTranslation,
                                                                           neighborIList, typeMolecule);
	}

    const double diffEnergy {newEnergyMolecule - oldEnergyMolecule};
    // Metropolis criterion
    const bool acceptMove{ metropolis(diffEnergy) };

    // If the move is accepted, then the energy, the position array and the displacement array can be updated.
    // If the m_calculatePressure is set to True, then the pressure is calculated.
    if (acceptMove)
    {
        generalUpdate(diffEnergy);
        /***
        if (m_calculatePressure)
        {
            const double newPressureParticle {pressureParticle(m_temp, indexTranslation, positionParticleTranslation, m_positionArray, neighborIList, m_typeArray, m_squareRc, m_lengthCube, m_halfLengthCube)};
            const double oldPressureParticle {pressureParticle(m_temp, indexTranslation, m_positionArray[indexTranslation], m_positionArray, neighborIList, m_typeArray, m_squareRc, m_lengthCube, m_halfLengthCube)};
            m_pressure += newPressureParticle - oldPressureParticle;
        }
        ***/
        m_systemMolecules.updateFlags(indexTranslation);
        m_systemMolecules.updatePositionI(indexTranslation, positionArrayTranslation);
        for (int j = 0; j < 3; j++)
        {
            // TO DO INTER DISPLACEMENT
            const int newIndexTranslation {indexTranslation + j };
            m_systemNeighbors.updateInterDisplacement(newIndexTranslation, randomVector);

        }
    }
    m_systemMolecules.reinitializeFlags();
}

/*******************************************************************************
 * This function returns a tentative new particle position.
 *
 * @param indexTranslation Translated particle's index.
 *        randomVector Random vector of size 3. The components of the vector are
 *                     taken from a uniform distribution U(-m_rBox, m_rBox).
 * @return Tentative new particle position.
 ******************************************************************************/
void MonteCarlo::mcTranslation()
{

    const int indexTranslation{randomIntGenerator(0, m_nParticles - 1)}; // randomly chosen particle
    const std::vector<double> randomVector(randomVectorDoubleGenerator(3, -m_rBox, m_rBox));
    //const std::vector<double>  randomVector = {0.02, 0, 0};
    std::vector<double> positionParticleTranslation { vectorTranslation(indexTranslation, randomVector) };
    const std::vector<int>& neighborIList { m_systemNeighbors.getNeighborIList(indexTranslation) };
    double oldEnergyParticle;
    double newEnergyParticle;

    oldEnergyParticle = m_systemMolecules.energyParticleMolecule(indexTranslation, neighborIList);
    newEnergyParticle = m_systemMolecules.energyParticleMolecule(indexTranslation, positionParticleTranslation,
                                                                 neighborIList);

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
        m_systemMolecules.updatePositionI(indexTranslation, positionParticleTranslation);
        m_systemMolecules.updateFlags(indexTranslation);
    }
    m_systemMolecules.reinitializeFlags();
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
    const std::vector<int>& neighborIList1 { m_systemNeighbors.getNeighborIList(indexSwap1) };
    const std::vector<int>& neighborIList2 { m_systemNeighbors.getNeighborIList(indexSwap2) };


    energyParticle1 = m_systemMolecules.energyParticleMolecule(indexSwap1, neighborIList1, indexSwap2);

    energyParticle2 = m_systemMolecules.energyParticleMolecule(indexSwap2, neighborIList2, indexSwap1);

    m_systemMolecules.swapParticleTypesIJ(indexSwap1, indexSwap2);

    energyParticleSwap1 = m_systemMolecules.energyParticleMolecule(indexSwap1, neighborIList1, indexSwap2);

    energyParticleSwap2 = m_systemMolecules.energyParticleMolecule(indexSwap2, neighborIList2, indexSwap1);


    const double diffEnergy{ energyParticleSwap1 + energyParticleSwap2 - energyParticle1 - energyParticle2 };
    // Metropolis criterion
    const bool acceptMove { metropolis( diffEnergy) };

    // If the move is accepted, then the energy, the position array and the displacement array can be updated.
    // If the m_calculatePressure is set to True, then the pressure is calculated.
    if ( acceptMove )
    {
        generalUpdate( diffEnergy );
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
        m_systemMolecules.swapParticleTypesIJ(indexSwap1, indexSwap2);
    }
}


std::vector<double> MonteCarlo::vectorTranslation(const int& indexTranslation, const std::vector<double>& randomVector)
{
    std::vector<double> positionParticleTranslation = m_systemMolecules.getPositionI(indexTranslation);
    positionParticleTranslation = vectorSum(positionParticleTranslation, randomVector);
    positionParticleTranslation = m_systemMolecules.periodicBC(positionParticleTranslation);
    return positionParticleTranslation;
}

void MonteCarlo::generalUpdate(double diffEnergy)
{
    m_energy += diffEnergy;
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