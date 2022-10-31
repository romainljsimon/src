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
};

struct inputVar
{
	double rc;
	double temp;
	double density;
	double rbox;
	double rskin;
	int neighUpdate;
	int numberIteration;
};


posRad readXYZ(std::string path);
void saveInXYZ(std::vector<std::vector<double>>& positionArray, std::vector<double> radiusArray, std::string path);
void saveEnergyTXT(double energy, std::string path);
inputVar readInput(std::string path);



#endif /* READSAVEFILE_H_ */
