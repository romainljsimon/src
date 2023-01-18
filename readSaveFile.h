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


posRad readXYZ(const std::string& path, const std::string& simulationMol);
void saveInXYZ(const std::vector<std::vector<double>>& positionArray, const std::vector<double>& radiusArray,
			   const std::vector<int>& moleculeType, const double& lengthCube,  const std::string& path);
void saveDoubleTXT(const double& number, const std::string& path);
void saveDisplacement(const std::vector<std::vector<double>>& dispMatrix, const std::string& path);
inputVar readInput(const std::string& path);
std::vector<std::vector<int>> readBondsTXT(const std::string& path);


#endif /* READSAVEFILE_H_ */
