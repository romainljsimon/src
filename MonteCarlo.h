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
#include "pressure.h"
#include "INPUT/Parameter.h"
#include "util.h"
#include "MOLECULES/Molecules.h"
#include "NEIGHBORS/Neighbors.h"
#include <cmath>


class MonteCarlo
{

private:
    Molecules m_systemMolecules;
    Neighbors m_systemNeighbors;
	double m_energy {};                                             // System's energy.
	double m_pressure {};                                           // System's pressure.
	const int m_nParticles {};                                            // System's number of particles.
    int m_nTrans {0};
	double m_acceptanceRateTrans { 0. };                                 // Monte Carlo translation acceptance rate.
    int m_nSwap {0};
    double m_acceptanceRateSwap { 0. };                                 // Monte Carlo swap acceptance rate.
    int m_nMolTrans {0};
    double m_acceptanceRateMolTrans { 0. };                                 // Monte Carlo swap acceptance rate.
    const int m_saveRate {};
	const bool m_calculatePressure {};                               // Boolean that decides if the pressure is calculated or not.
    const bool m_swap {};
    const double m_pSwap {};
	const std::string m_simulationMol {};                       			// Type of system: can be either "polymer" or "atomic".
    const bool m_molTranslation {};
    const double m_pMolTranslation {};
    const double m_rBoxMolTrans {};
	const double m_temp {};                                     	// Temperature.
	const double m_rBox{};                               			// Length of the translation box.
	const int m_saveUpdate {};                           			// save xyz update frequency.
	const std::string m_neighMethod {};                   			// Neighbor list method: "verlet" for verlet neighbor list. Any other value: no neighbor list.
	const int m_timeSteps {};                             			// Number of time steps.
    const std::string m_folderPath {};                    			// Path to where the algorithm was launched.


public:
    // Monte Carlo constructor
    MonteCarlo (param::Parameter param, const Molecules& systemMolecules, Neighbors  systemNeighbors,
                std::string folderPath)

        : m_systemMolecules (systemMolecules)
            , m_systemNeighbors(std::move(systemNeighbors))
            , m_nParticles(systemMolecules.getNParticles())
            , m_calculatePressure(param.get_bool("calcPressure", false))
            , m_swap (param.get_bool("swap", false))
            , m_pSwap (param.get_double("pSwap", 0.2))
            , m_molTranslation ( param.get_bool("molTranslation", false))
            , m_pMolTranslation ( param.get_double("pMolTranslation", 0.1))
            , m_rBoxMolTrans ( param.get_double("rBoxMolTranslation", 0.05))
            , m_temp { param.get_double( "temp") }
            , m_rBox { param.get_double( "rBox") }
            , m_saveUpdate { param.get_int( "waitingTime") }
            , m_timeSteps { param.get_int( "timeSteps") }
            , m_saveRate { param.get_int("saveRate", 1000)}
            , m_folderPath (std::move( folderPath ))

    {
        m_energy = m_systemMolecules.energySystemMolecule( m_systemNeighbors );
    }

	void mcTotal();
	int mcMove();
    void mcTranslation();
    void generalUpdate(double diff_energy);
    void mcSwap();
    [[nodiscard]] bool metropolis(double diff_energy) const;

    template<typename InputIt>
    std::vector<double> vectorTranslation(const int& indexTranslation, InputIt randomVectorIt)
    {
        auto posItBeginTranslation = m_systemMolecules.getPosItBeginI(indexTranslation);
        const int& nDims {m_systemMolecules.getNDims()};
        std::vector<double> positionTranslation (nDims);
        std::transform(posItBeginTranslation, posItBeginTranslation + nDims,
                       randomVectorIt, positionTranslation.begin(), std::plus<>());
        m_systemMolecules.periodicBC(positionTranslation.begin());
        return positionTranslation;
    }


    void mcMoleculeTranslation();
};


#endif /* MONTECARLO_H_ */
