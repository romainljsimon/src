/*
 * MonteCarlo.h
 *
 *  Created on: 27 oct. 2022
 *      Author: Romain Simon
 */

#ifndef MONTECARLO_H_
#define MONTECARLO_H_

#include <utility>
#include "energy.h"
#include "pressure.h"
#include "input/Parameter.h"

class MonteCarlo
{

private:
	std::vector<std::vector<int>> m_neighborList {};                // Neighbor list.
	int m_errors { 0 };                                             // Errors of the neighbor list.
	double m_energy {};                                             // System's energy.
	double m_pressure {};                                           // System's pressure.
	const int m_nParticles {};                                            // System's number of particles.
	std::vector<std::vector<double>> m_totalDisplacementMatrix {};  // Total displacement Matrix.
	std::vector<std::vector<double>> m_interDisplacementMatrix {};  // Inter neighbor list update displacement matrix.
	std::vector<std::vector<double>> m_stepDisplacementMatrix {};   // Step displacement matrix.
	double m_acceptanceRate { 0. };                                 // Monte Carlo acceptance rate.
    double m_acceptanceRateSwap { 0. };                                 // Monte Carlo acceptance rate.
    double m_updateRate { -1. };                                     // Monte Carlo neighbor list update rate.
    const int m_saveRate {};
	const bool m_calculatePressure {};                               // Boolean that decides if the pressure is calculated or not.
    const bool m_swap {};
    const double m_pSwap {};
	const std::string m_simulationMol {};                       			// Type of system: can be either "polymer" or "atomic".
    const double m_lengthCube {};                                   // Length of the simulation box.
    const double m_rC {};
	const double m_squareRc {};                           			// Cut off radius squared.

	const double m_temp {};                                     	// Temperature.
	const double m_rBox{};                               			// Length of the translation box.
    const double m_rSkin {};
	const double m_squareRSkin {};                        			// Skin radius squared.
	const int m_saveUpdate {};                           			// save xyz update frequency.
	const std::string m_neighMethod {};                   			// Neighbor list method: "verlet" for verlet neighbor list. Any other value: no neighbor list.
	const int m_timeSteps {};                             			// Number of time steps.
	const double m_squareR0 {};                           			// If the simulation is polymeric: max length of a bond.
	const double m_feneK {};                              			// If the simulation is polymeric: stiffness of a bond.
	const double m_squareRDiff {};                        			// squared difference of skin - cut off.
    std::vector<std::vector<double>> m_positionArray {};
    const std::vector<std::vector<int>> m_bondsMatrix {};           // If the simulation is polymeric: matrix of bonded nearest-neighbors.// Particles positions array of size (N, 3).
    std::vector<double> m_diameterArray {};                 		// Particles diameters array of size (N, 1).
    const std::vector<int> m_moleculeType {};                       // Particles type array of size (N, 1).
    const std::string m_folderPath {};                    			// Path to where the algorithm was launched.


public:
    // Polymer constructor
    MonteCarlo (param::Parameter param, std::vector<std::vector<double>> positionArray,
                std::vector<double> diameterArray, std::vector<int> moleculeType,
                std::vector<std::vector<int>> bondsMatrix, std::string folderPath)

            : m_simulationMol (param.get_string("simType"))
            , m_calculatePressure(param.get_bool("calcPressure", false))
            , m_swap (param.get_bool("swap", false))
            , m_pSwap (param.get_double("pSwap", 0.2))
            , m_rC { param.get_double( "rc") }
            , m_squareRc { pow (param.get_double( "rc"), 2 ) }
            , m_lengthCube { pow (static_cast<double>( positionArray.size()) / param.get_double("density"), 1. / 3.)}
            , m_temp { param.get_double( "temp") }
            , m_rBox { param.get_double( "rBox") }
            , m_rSkin { param.get_double( "rSkin") }
            , m_squareRSkin {pow (param.get_double( "rSkin"), 2 ) }
            , m_saveUpdate { param.get_int( "waitingTime") }
            , m_neighMethod (param.get_string( "neighMethod", "verlet"))
            , m_timeSteps { param.get_int( "timeSteps") }
            , m_saveRate { param.get_int("saveRate", 50)}
            , m_squareR0 { pow (param.get_double( "r0", 1.5), 2 ) }
            , m_feneK { param.get_double( "feneK", 30.)}
            , m_squareRDiff{pow ((param.get_double( "rSkin") - param.get_double(
                     "rc")) / 2, 2)}
            , m_nParticles {static_cast<int>( positionArray.size() )}
            , m_positionArray (std::move( positionArray ))
            , m_diameterArray (std::move( diameterArray ))
            , m_moleculeType (std::move( moleculeType ))
            , m_bondsMatrix (std::move(bondsMatrix))
            , m_folderPath (std::move( folderPath ))


    {
        m_interDisplacementMatrix.resize(m_nParticles, std::vector<double>(3, 0));
        m_totalDisplacementMatrix.resize(m_nParticles, std::vector<double>(3, 0));
        m_stepDisplacementMatrix.resize(m_nParticles, std::vector<double>(3, 0));
        m_neighborList.resize(m_nParticles);

        createNeighborList();

        m_energy = energySystemPolymer(m_positionArray, m_diameterArray, m_bondsMatrix,
                                           m_neighborList, m_squareRc, m_lengthCube, m_squareR0, m_feneK);
        if (m_calculatePressure)
        {
            m_pressure = pressureSystem(m_temp, m_positionArray, m_diameterArray, m_squareRc, m_lengthCube); //+ corrPressure;
        }

    }
    // Atomic constructor
    MonteCarlo (param::Parameter param, std::vector<std::vector<double>> positionArray,
                std::vector<double> diameterArray, std::vector<int> moleculeType, std::string folderPath)


            : m_simulationMol (param.get_string( "simType"))
            , m_calculatePressure(param.get_bool("calcPressure", false))
            , m_swap (param.get_bool("swap", false))
            , m_pSwap (param.get_double("pSwap", 0.2))
            , m_rC { param.get_double( "rc") }
            , m_squareRc { pow (param.get_double( "rc"), 2 ) }
            , m_lengthCube { pow (static_cast<double>(positionArray.size()) / param.get_double(
                     "density"), 1. / 3)}
            , m_temp { param.get_double( "temp") }
            , m_rBox { param.get_double( "rBox") }
            , m_rSkin { param.get_double( "rSkin") }
            , m_squareRSkin {pow (param.get_double( "rSkin"), 2 ) }
            , m_saveUpdate { param.get_int( "waitingTime") }
            , m_neighMethod (param.get_string( "neighMethod", "verlet"))
            , m_timeSteps { param.get_int( "timeSteps") }
            , m_saveRate { param.get_int("saveRate", 50)}
            , m_squareR0 { pow (param.get_double( "r0", 1.5), 2 ) }
            , m_feneK { param.get_double( "feneK", 30.)}
            , m_nParticles {static_cast<int>( positionArray.size() )}
            , m_positionArray (std::move( positionArray ))
            , m_diameterArray (std::move( diameterArray ))
            , m_moleculeType (std::move( moleculeType ))
            , m_folderPath (std::move( folderPath ))


    {
        m_interDisplacementMatrix.resize(m_nParticles, std::vector<double>(3, 0));
        m_totalDisplacementMatrix.resize(m_nParticles, std::vector<double>(3, 0));
        m_stepDisplacementMatrix.resize(m_nParticles, std::vector<double>(3, 0));
        m_neighborList.resize(m_nParticles);

        createNeighborList();

        m_energy = energySystem( m_positionArray, m_diameterArray,
                                 m_neighborList, m_squareRc, m_lengthCube);
        if (m_calculatePressure)
        {
            m_pressure = pressureSystem(m_temp, m_positionArray, m_diameterArray,
                                        m_squareRc, m_lengthCube); //+ corrPressure;
        }

    }

	void mcTotal();
	void createNeighborList();
	void mcMove();
	void checkStepDisplacement();
    std::vector<int> findNeighborIList(int indexTranslation);
    void mcTranslation();
    void generalUpdate(double diff_energy);
    void mcSwap();
    [[nodiscard]] bool metropolis(double diff_energy) const;
};


#endif /* MONTECARLO_H_ */
