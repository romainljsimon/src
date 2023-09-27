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
#include "Random_mt.h"
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
    const std::string energyFilePath{"./outE.txt"};
    //std::string pressureFilePath {"./outP.txt"};
    // std::vector<double> radiusArray (divideVectorByScalar(m_typeArray, 2));

	m_systemMolecules.saveInXYZ(preName + std::to_string(0) + extname );
	saveDoubleTXT(m_energy / m_nParticles, energyFilePath);
	// saveDisplacement(m_totalDisplacementMatrix, preNameDisp + std::to_string(0) + extnameDisp);

	const std::vector<int> saveTimeStepArray ( createSaveTime(m_timeSteps, m_saveUpdate, 1.1));

	int save_index { 0 };
    double swapRate {0.};
    double transRate {0.};
    double molTransRate {0.};
	for (int i = 0; i < m_timeSteps; i++) //Iteration over m_timeSteps
	{
        int j { 0 };
        while (j < m_nParticles) // N Monte Carlo moves are tried in one time step.
		{
            j += mcMove();
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
        swapRate += static_cast<double>(m_nSwap) / m_nParticles;
        transRate += static_cast<double>(m_nTrans) / m_nParticles;
        molTransRate += static_cast<double>(m_nMolTrans) / m_nParticles;
        m_nSwap = 0;
        m_nTrans = 0;
        m_nMolTrans = 0;

        /***
		if (m_calculatePressure) // Saves pressure if m_calculatePressure is True
		{
			saveDoubleTXT( m_pressure, pressureFilePath);
		}
        ***/

	}

    m_systemMolecules.saveInXYZ( preName + std::to_string(m_timeSteps) + extname );
    const double totalRate { transRate + swapRate + molTransRate};
    double totalAcceptanceRate = (totalRate != 0) ?
            (m_acceptanceRateTrans + m_acceptanceRateSwap + m_acceptanceRateMolTrans) / totalRate : 0.;
	m_acceptanceRateTrans = (transRate != 0.) ? m_acceptanceRateTrans / transRate : 0.;
    m_acceptanceRateSwap = (swapRate != 0.) ? m_acceptanceRateSwap / swapRate : 0.;
    m_acceptanceRateMolTrans = (molTransRate != 0.) ? m_acceptanceRateMolTrans / molTransRate : 0.;


    constexpr std::string_view translationString{"Translation MC move acceptance rate: "};
    constexpr std::string_view swapString { "Swap MC move acceptance rate: " };
    constexpr std::string_view molTranslationString { "Molecule translation MC move acceptance rate: " };
    constexpr std::string_view totalString { "Total MC move acceptance rate: " };

    std::cout << translationString << m_acceptanceRateTrans << "\n";
    std::cout << swapString << m_acceptanceRateSwap << "\n";
    std::cout << molTranslationString << m_acceptanceRateMolTrans << "\n";
    std::cout << totalString << totalAcceptanceRate << "\n";


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
int MonteCarlo::mcMove()
{
    const double randomDouble { Random::doubleGenerator(0., 1.) } ;
    const bool swapped {randomDouble < m_pSwap};
    const bool molTranslation {randomDouble > (1 - m_pMolTranslation)};
    int step {0};

    if ( swapped && m_swap)
    {

        ++m_nSwap;
        ++step;
        mcSwap();
    }
    else if ( molTranslation && m_molTranslation )
    {
        ++m_nMolTrans;
        step += 3;
        mcMoleculeTranslation();
    }
    else
    {

        ++m_nTrans;
        ++step;
        mcTranslation();
    }
    return step;
}

/*******************************************************************************
 * This function implements a Monte Carlo move: translation of a the whole
 * polymer, calculation of the energy of the new system and then acceptation or not of
 * the move according to the Metropolis criterion. It is called ??? times in one
 * time step.
 ******************************************************************************/

void MonteCarlo::mcMoleculeTranslation()
{
	constexpr int lenMolecule {3};
    const int typeMolecule { Random::intGenerator(0, (m_nParticles - 1) / lenMolecule) }; // randomly chosen Molecule
    const int indexTranslation { m_systemMolecules.getNDims() * typeMolecule};
	std::vector<double> randomVector ( Random::vectorDoubleGenerator(3, -m_rBoxMolTrans, m_rBoxMolTrans) );
	double oldEnergyMolecule {0};
	double newEnergyMolecule {0};
    std::vector<double> positionArrayTranslation;

	for (int j = 0; j < lenMolecule; j++)
	{

        const int newIndexTranslation {indexTranslation + j };
        std::vector<double> posTranslation{ vectorTranslation(newIndexTranslation,
                                                            randomVector.begin())};

        positionArrayTranslation.insert( positionArrayTranslation.end(),
                                         posTranslation.begin(),
                                         posTranslation.begin() + m_systemMolecules.getNDims());

        const auto& neighItBegin { m_systemNeighbors.getNeighItBeginI(newIndexTranslation) };
        const int& lenNeigh {m_systemNeighbors.getLenIndexBegin(newIndexTranslation)};

        oldEnergyMolecule += m_systemMolecules.energyPairParticleExtraMolecule( newIndexTranslation,
                                                                                neighItBegin, lenNeigh,
                                                                                typeMolecule);
		newEnergyMolecule += m_systemMolecules.energyPairParticleExtraMolecule( newIndexTranslation,
                                                                                posTranslation.begin(),
                                                                                neighItBegin, lenNeigh,
                                                                                typeMolecule);
	}

    const double diffEnergy {newEnergyMolecule - oldEnergyMolecule};
    // Metropolis criterion
    const bool acceptMove{ metropolis(diffEnergy) };

    // If the move is accepted, then the energy, the position array and the displacement array can be updated.
    // If the m_calculatePressure is set to True, then the pressure is calculated.
    if (acceptMove)
    {
        generalUpdate(diffEnergy);
        m_acceptanceRateMolTrans += 1. / m_nParticles; // increment of the acceptance rate.


        //if (m_calculatePressure)
        //{
        //    const double newPressureParticle {pressureParticle(m_temp, indexTranslation, positionParticleTranslation, m_positionArray, neighborIList, m_typeArray, m_squareRc, m_lengthCube, m_halfLengthCube)};
        //    const double oldPressureParticle {pressureParticle(m_temp, indexTranslation, m_positionArray[indexTranslation], m_positionArray, neighborIList, m_typeArray, m_squareRc, m_lengthCube, m_halfLengthCube)};
        //    m_pressure += newPressureParticle - oldPressureParticle;
        //}

        m_systemMolecules.updateFlags(indexTranslation);
        m_systemMolecules.updatePositionI(indexTranslation, positionArrayTranslation.begin(), lenMolecule*m_systemMolecules.getNDims());
        for (int j = 0; j < lenMolecule; j++)
        {
            // TO DO INTER DISPLACEMENT
            const int newIndexTranslation {indexTranslation + j };
            m_systemNeighbors.updateInterDisplacement(newIndexTranslation,
                                                      randomVector.begin());

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

    const int indexTranslation{Random::intGenerator(0, m_nParticles - 1)}; // randomly chosen particle
    const std::vector<double>& randomVector(Random::vectorDoubleGenerator(3, -m_rBox, m_rBox));
    const std::vector<double>& positionTranslation { vectorTranslation(indexTranslation, randomVector.begin()) };

    const auto& neighItBegin { m_systemNeighbors.getNeighItBeginI(indexTranslation) };
    const int& lenNeigh {m_systemNeighbors.getLenIndexBegin(indexTranslation)};

    double oldEnergyParticle { m_systemMolecules.energyParticleMolecule(indexTranslation, neighItBegin,
                                                                        lenNeigh)};
    double newEnergyParticle { m_systemMolecules.energyParticleMolecule(indexTranslation,
                                                                        positionTranslation.begin(),
                                                                        neighItBegin, lenNeigh)};

    const double diff_energy {newEnergyParticle - oldEnergyParticle};

    // Metropolis criterion
    const bool acceptMove{ metropolis(diff_energy) };

    // If the move is accepted, then the energy, the position array and the displacement array can be updated.
    // If the m_calculatePressure is set to True, then the pressure is calculated.
    if (acceptMove)
    {
        generalUpdate(diff_energy);
        m_acceptanceRateTrans += 1. / m_nParticles; // increment of the acceptance rate.

        /***
        if (m_calculatePressure)
        {
            const double newPressureParticle {pressureParticle(m_temp, indexTranslation, positionParticleTranslation, m_positionArray, neighborIList, m_typeArray, m_squareRc, m_lengthCube, m_halfLengthCube)};
            const double oldPressureParticle {pressureParticle(m_temp, indexTranslation, m_positionArray[indexTranslation], m_positionArray, neighborIList, m_typeArray, m_squareRc, m_lengthCube, m_halfLengthCube)};
            m_pressure += newPressureParticle - oldPressureParticle;
        }
        ***/
        m_systemNeighbors.updateInterDisplacement(indexTranslation, randomVector.begin());
        m_systemMolecules.updatePositionI(indexTranslation, positionTranslation.begin());
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
    int indexMolecule { Random::intGenerator(0, m_nParticles - 1) };
    // int indexSwap2 { randomIntGenerator(0, m_nParticles - 1) };
    // !!!! THIS NEXT PART IS ONLY BECAUSE WE STUDY TRI-MERS VERY SPECIFIC!!!
    indexMolecule -= indexMolecule % 3;
    int indexSwap1 { Random::intGenerator(0, 2)};
    int indexSwap2 { Random::intGenerator(0, 2)};

    while (indexSwap1 == indexSwap2)
    {
        indexSwap2 = Random::intGenerator(0, 2);
    }


    indexSwap1 += indexMolecule;
    indexSwap2 += indexMolecule;

    const auto& neighItBegin1 { m_systemNeighbors.getNeighItBeginI(indexSwap1) };
    const int& lenNeigh1 {m_systemNeighbors.getLenIndexBegin(indexSwap1)};
    const auto& neighItBegin2 { m_systemNeighbors.getNeighItBeginI(indexSwap2) };
    const int& lenNeigh2 {m_systemNeighbors.getLenIndexBegin(indexSwap2)};


    double diffEnergySwap1{ m_systemMolecules.energyParticleMoleculeSwap(indexSwap1, neighItBegin1,
                                                                   lenNeigh1, indexSwap2)};


    double diffEnergySwap2 { m_systemMolecules.energyParticleMoleculeSwap(indexSwap2, neighItBegin2,
                                                               lenNeigh2, indexSwap1)};




    const double diffEnergy{ diffEnergySwap1 + diffEnergySwap2 };
    // Metropolis criterion
    const bool acceptMove { metropolis( diffEnergy) };

    // If the move is accepted, then the energy, the position array and the displacement array can be updated.
    // If the m_calculatePressure is set to True, then the pressure is calculated.
    if ( acceptMove )
    {
        generalUpdate( diffEnergy );
        m_acceptanceRateSwap += 1. / m_nParticles;
        m_systemMolecules.swapParticleTypesIJ(indexSwap1, indexSwap2);

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
}


void MonteCarlo::generalUpdate(double diffEnergy)
{
    m_energy += diffEnergy;

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
		const double randomDouble { Random::doubleGenerator(0., 1.) } ;
		const double threshold { exp(( - diff_energy ) / m_temp) }; // we consider k=1
		return threshold > randomDouble;
	}
}