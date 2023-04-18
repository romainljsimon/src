/*
 * readSaveFile.h
 *
 *  Created on: 7 oct. 2022
 *      Author: Romain Simon
 */

#ifndef READSAVEFILE_H_
#define READSAVEFILE_H_

struct posRad
{
	std::vector<std::vector<double>> posMatrix;
	std::vector<double> radVector;
	std::vector<int> moleculeType;
};



posRad readXYZ(const std::string& path);
void saveInXYZ(const std::vector<std::vector<double>>& positionArray, const std::vector<double>& radiusArray,
			   const std::vector<int>& moleculeType, const double& lengthCube,  const std::string& path);
void saveDoubleTXT(const double& number, const std::string& path);
void saveDisplacement(const std::vector<std::vector<double>>& dispMatrix, const std::string& path);
std::vector<std::vector<int>> readBondsTXT(const std::string& path);


#endif /* READSAVEFILE_H_ */
