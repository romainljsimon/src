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
	std::vector<std::vector<int>> m_neighborList {};
	int m_errors { 0 };
	std::vector<std::vector<std::vector<std::vector<int>>>> m_cellList {};
	double m_energy {};
	double m_pressure {};
	int m_nParticles {};
	std::vector<std::vector <int>> m_particleIndexCell {};
	std::vector<std::vector<int>> m_bondsMatrix {};
	std::vector<std::vector<double>> m_totalDisplacementMatrix {};
	std::vector<std::vector<double>> m_interDisplacementMatrix {};
	std::vector<std::vector<double>> m_stepDisplacementMatrix {};
	double m_acceptanceRate { 0 };
	double m_updateRate { -1 };
	bool m_calculatePressure {false};

	std::string m_simulationMol {};
	std::vector<std::vector<double>> m_positionArray {};
	std::vector<double> m_radiusArray {};
	std::vector<int> m_moleculeType {};
	const double m_squareRc {};
	const double m_lengthCube {};
	double m_temp {};
	const double m_rbox {};
	const double m_squareRskin {};
	const int m_neighUpdate {};
	const std::string m_folderPath {};
	const std::string m_neighMethod {};
	const int m_timeSteps {};
	const double m_squareR0 {};
	const double m_feneK {};
	const double m_squareRdiff {};

public:
	MonteCarlo( std::string simulationMol, std::vector<std::vector<double>> positionArray,
				std::vector<double> radiusArray, std::vector<int> moleculeType,
				const double rc, const double lengthCube, const double temp, const double rbox,
				const double rskin, const int neighUpdate, const std::string folderPath,
				const std::string neighMethod, const int timeSteps, const double r0,
				const double feneK);

	void mcTotal();
	void createNeighborList();
	void createCellAndIndexList();
	void mcMove();
	std::vector<double> mcTranslation(int indexTranslation, std::vector<double> randomVector);
	bool metropolis(double newEnergy, double energy);
	void checkStepDisplacement();

};


#endif /* MONTECARLO_H_ */
