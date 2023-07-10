/*
 * MonteCarlo.h
 *
 *  Created on: 27 oct. 2022
 *      Author: Romain Simon
 */

#ifndef MONTECARLO_H_
#define MONTECARLO_H_

#include <utility>
#include <iterator>
#include <fstream>
#include <vector>
#include "energy.h"
#include "pressure.h"
#include "INPUT/Parameter.h"
#include "util.h"
#include "POTENTIALS/PairPotentials.h"
#include "POTENTIALS/BondPotentials.h"
#include "PARTICLES/Particles.h"
#include "NEIGHBORS/Neighbors.h"


class MonteCarlo
{

private:
    const BondPotentials m_bondPotentials;
    const PairPotentials m_pairPotentials;
    Particles m_systemParticles;
    Neighbors m_systemNeighbors;
	int m_errors { 0 };                                             // Errors of the neighbor list.
	double m_energy {};                                             // System's energy.
	double m_pressure {};                                           // System's pressure.
	const int m_nParticles {};                                            // System's number of particles.
	double m_acceptanceRate { 0. };                                 // Monte Carlo acceptance rate.
    double m_acceptanceRateSwap { 0. };                                 // Monte Carlo acceptance rate.
    double m_updateRate { -1. };                                     // Monte Carlo neighbor list update rate.
    const int m_saveRate {};
	const bool m_calculatePressure {};                               // Boolean that decides if the pressure is calculated or not.
    const bool m_swap {};
    const double m_pSwap {};
	const std::string m_simulationMol {};                       			// Type of system: can be either "polymer" or "atomic".

	const double m_temp {};                                     	// Temperature.
	const double m_rBox{};                               			// Length of the translation box.
	const int m_saveUpdate {};                           			// save xyz update frequency.
	const std::string m_neighMethod {};                   			// Neighbor list method: "verlet" for verlet neighbor list. Any other value: no neighbor list.
	const int m_timeSteps {};                             			// Number of time steps.
    const std::string m_folderPath {};                    			// Path to where the algorithm was launched.


public:
    // Polymer constructor
    MonteCarlo (param::Parameter param, const Particles& systemParticles, const PairPotentials& systemPairPotentials,
                const BondPotentials& systemBondPotentials, const Neighbors& systemNeighbors, std::string folderPath)

            : m_simulationMol (param.get_string("simType"))
            , m_bondPotentials (systemBondPotentials)
            , m_pairPotentials (systemPairPotentials)
            , m_systemParticles (systemParticles)
            , m_systemNeighbors(systemNeighbors)
            , m_nParticles(systemParticles.getNParticles())
            , m_calculatePressure(param.get_bool("calcPressure", false))
            , m_swap (param.get_bool("swap", false))
            , m_pSwap (param.get_double("pSwap", 0.2))
            , m_temp { param.get_double( "temp") }
            , m_rBox { param.get_double( "rBox") }
            , m_saveUpdate { param.get_int( "waitingTime") }
            , m_neighMethod (param.get_string( "neighMethod", "verlet"))
            , m_timeSteps { param.get_int( "timeSteps") }
            , m_saveRate { param.get_int("saveRate", 50)}
            , m_folderPath (std::move( folderPath ))
            , m_energy { energySystemPolymer( systemParticles, systemNeighbors, systemPairPotentials, systemBondPotentials) }


    {
        // Potential Initializations


        //if (m_calculatePressure)
        //{
        //    m_pressure = pressureSystem(m_temp, m_positionArray, m_typeArray, m_squareRc, m_lengthCube, m_halfLengthCube); //+ corrPressure;
        //}

    }
    // Atomic constructor
    MonteCarlo (param::Parameter param, const Particles& systemParticles,
                const PairPotentials& systemPairPotentials,  const Neighbors& systemNeighbors, std::string folderPath)

            : m_simulationMol (param.get_string("simType"))
            , m_pairPotentials(systemPairPotentials)
            , m_systemParticles (systemParticles)
            , m_calculatePressure(param.get_bool("calcPressure", false))
            , m_nParticles(systemParticles.getNParticles())
            , m_swap (param.get_bool("swap", false))
            , m_pSwap (param.get_double("pSwap", 0.2))
            , m_systemNeighbors(systemNeighbors)
            , m_temp { param.get_double( "temp") }
            , m_rBox { param.get_double( "rBox") }
            , m_saveUpdate { param.get_int( "waitingTime") }
            , m_neighMethod (param.get_string( "neighMethod", "verlet"))
            , m_timeSteps { param.get_int( "timeSteps") }
            , m_saveRate { param.get_int("saveRate", 50)}
            , m_folderPath (std::move( folderPath ))
            , m_energy { energySystem( systemParticles, systemNeighbors, systemPairPotentials) }
    {



        //if (m_calculatePressure)
        //{
        //    m_pressure = pressureSystem(m_temp, m_positionArray, m_typeArray,
         //                               m_squareRc, m_lengthCube, m_halfLengthCube); //+ corrPressure;
        //}

    }

	void mcTotal();

	void mcMove();
    void mcTranslation();
    void generalUpdate(double diff_energy);
    void mcSwap();
    [[nodiscard]] bool metropolis(double diff_energy) const;




};


#endif /* MONTECARLO_H_ */
