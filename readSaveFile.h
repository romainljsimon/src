/*
 * readSaveFile.h
 *
 *  Created on: 7 oct. 2022
 *      Author: Romain Simon
 */

#ifndef READSAVEFILE_H_
#define READSAVEFILE_H_


#include <vector>

void saveInXYZ(const std::vector<std::vector<double>>& positionArray, const std::vector<double>& radiusArray,
               const std::vector<int>& moleculeType, const double& lengthCube,  const std::string& path);
void saveDoubleTXT(const double& number, const std::string& path);
void saveDisplacement(const std::vector<std::vector<double>>& dispMatrix, const std::string& path);
std::vector<std::vector<int>> readBondsTXT(const std::string& path);
void saveDoubleIntTXT(const double& number1, const int& number2, const std::string& path);


#endif /* READSAVEFILE_H_ */
