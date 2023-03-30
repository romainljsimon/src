/*
 * MonteCarlo.h
 *
 *  Created on: 27 oct. 2022
 *      Author: rsimon
 */

#ifndef MONTECARLO_H_
#define MONTECARLO_H_



class MonteCarlo
{

private:
	std::vector<std::vector<int>> m_neighborList {};                // Neighbor list.
	int m_errors { 0 };                                             // Errors of the neighbor list.
	double m_energy {};                                             // System's energy.
	double m_pressure {};                                           // System's pressure.
	int m_nParticles {};                                            // System's number of particles.
	std::vector<std::vector<int>> m_bondsMatrix {};                 // If the simulation is polymeric: matrix of bonded nearest-neighbors.
	std::vector<std::vector<double>> m_totalDisplacementMatrix {};  // Total displacement Matrix.
	std::vector<std::vector<double>> m_interDisplacementMatrix {};  // Inter neighbor list update displacement matrix.
	std::vector<std::vector<double>> m_stepDisplacementMatrix {};   // Step displacement matrix.
	double m_acceptanceRate { 0. };                                 // Monte Carlo acceptance rate.
	double m_updateRate { -1. };                                     // Monte Carle neighbor list update rate.
	bool m_calculatePressure {false};                               // Boolean that decides if the pressure is calculated or not.

	std::string m_simulationMol {};                       			// Type of system: can be either "polymer" or "atomic".
	std::vector<std::vector<double>> m_positionArray {};  			// Particles positions array of size (N, 3).
	std::vector<double> m_diameterArray {};                 		// Particles diameters array of size (N, 1).
	std::vector<int> m_moleculeType {};                   			// Particles type array of size (N, 1).
	const double m_squareRc {};                           			// Cut off radius squared.
	const double m_lengthCube {};                         			// Length of the simulation box.
	double m_temp {};                                     			// Temperature.
	const double m_rbox {};                               			// Length of the translation box.
	const double m_squareRskin {};                        			// Skin radius squared.
	const int m_saveUpdate {};                           			// save xyz update frequency.
	const std::string m_folderPath {};                    			// Path to where the algorithm was launched.
	const std::string m_neighMethod {};                   			// Neighbor list method: "verlet" for verlet neighbor list. Any other value: no neighbor list.
	const int m_timeSteps {};                             			// Number of time steps.
	const double m_squareR0 {};                           			// If the simulation is polymeric: max length of a bond.
	const double m_feneK {};                              			// If the simulation is polymeric: stiffness of a bond.
	const double m_squareRdiff {};                        			// squared difference of skin - cut off.

public:
	MonteCarlo( std::string simulationMol, std::vector<std::vector<double>> positionArray,
				std::vector<double> diameterArray, std::vector<int> moleculeType,
				double rc, double lengthCube, double temp, double rbox,
				double rskin, int saveUpdate, std::string  folderPath,
				std::string  neighMethod, int timeSteps, double r0,
				double feneK);

	void mcTotal();
	void createNeighborList();
	// void mcAllMove();
	void mcMove();
	std::vector<double> mcTranslation(int indexTranslation, const std::vector<double>& randomVector);
	bool metropolis(double newEnergy, double energy);
	void checkStepDisplacement();

};


#endif /* MONTECARLO_H_ */
