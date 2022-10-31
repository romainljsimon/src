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
	int m_nParticles {};
	std::vector<std::vector <int>> m_particleIndexCell {};
	bool m_acceptMove {};
	int m_indexTranslation {};
	std::vector<double> m_positionParticleTranslation;

	std::vector<std::vector<double>> m_positionArray {};
	std::vector<double> m_radiusArray {};
	double m_rc {};
	double m_lengthCube {};
	double m_temp {};
	double m_rbox {};
	double m_rskin {};
	int m_neighUpdate {};
	int m_numberIteration {};
	std::string m_folderPath {};
	std::string m_neighMethod { };

public:
	MonteCarlo( std::vector<std::vector<double>> positionArray, std::vector<double> radiusArray,
				double rc, double lengthCube, double temp, double rbox, double rskin, int neighUpdate,
				std::string folderPath, std::string neighMethod, int numberIteration );

	void mcTotal();
	void createNeighborList();
	void createCellAndIndexList();
	void mcMove();
	void mcTranslation();
	void metropolis(double newEnergy, double energy);
};


#endif /* MONTECARLO_H_ */
