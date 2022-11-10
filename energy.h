/*
 * energy.h
 *
 *  Created on: 3 oct. 2022
 *      Author: romainsimon
 */

#ifndef ENERGY_H_
#define ENERGY_H_

struct cellAndIndex
{
	std::vector<std::vector<std::vector<std::vector<int>>>> cellList;
	std::vector<std::vector <int>> particleIndexCell;
};

double ljPotential(double distance, double sigmaA, double sigmaB, double squareRc);
double energySystem(std::vector<std::vector<double>> positionArray, std::vector<double> radiusArray, double squareRc, double lengthCube);
double squareDistancePair(std::vector<double> positionA, std::vector<double> positionB, double lengthCube);
double energyParticle(int indexParticle, std::vector<double> positionParticle, std::vector<std::vector<double>> positionArray, std::vector<int> neighborIList, std::vector<double> radiusArray, double squareRc, double lengthCube);
double energyParticlePolymer (int indexParticle, std::vector<double> positionParticle,
								std::vector<std::vector<double>> positionArray, std::vector<int> neighborIList,
								std::vector<double> radiusArray, std::vector<int> bondsI,
								double squareRc, double lengthCube, double squareR0, double feneK);

double energySystemPolymer(std::vector<std::vector<double>> positionArray, std::vector<double> radiusArray,
						   std::vector<std::vector<int>> bondsMatrix, double squareRc, double lengthCube,
						   double squareR0, double feneK);
std::vector<std::vector<int>> createNeighborList(std::vector<std::vector<double>> positionArray, double skin, double lengthCube);
std::vector<std::vector<std::vector<std::vector<int>>>> createCellList(std::vector<std::vector<double>> positionArray, double skin, double lengthCube);
double pressureSystem(double temp, std::vector<std::vector<double>> positionArray, std::vector<double> radiusArray, double lengthCube);





#endif /* ENERGY_H_ */
