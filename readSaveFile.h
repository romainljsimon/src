/*
 * readSaveFile.h
 *
 *  Created on: 7 oct. 2022
 *      Author: romainsimon
 */

#ifndef READSAVEFILE_H_
#define READSAVEFILE_H_

struct posRad
{
	std::vector<std::vector<double>> posMatrix;
	std::vector<double> radVector;
	std::vector<int> moleculeType;
};

struct inputVar
{
	std::string simulationMol;
	std::string neighborMethod;
	double rc;
	double temp;
	double density;
	double rbox;
	double rskin;
	int neighUpdate;
	int timeSteps;
	double r0;
	double feneK;
};


posRad readXYZ(std::string path, std::string simulationMol);
void saveInXYZ(std::vector<std::vector<double>>& positionArray, std::vector<double> radiusArray,
			   std::vector<int> moleculeType,double lengthCube,  std::string path);
void saveEnergyTXT(double energy, std::string path);
void saveDisplacement(std::vector<std::vector<double>> dispMatrix, std::string path);

inputVar readInput(std::string path);
std::vector<std::vector<int>> readBondsTXT(std::string path);



#endif /* READSAVEFILE_H_ */
