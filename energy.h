/*
 * energy.h
 *
 *  Created on: 3 oct. 2022
 *      Author: romainsimon
 */

#ifndef ENERGY_H_
#define ENERGY_H_

double ljPotential(double distance, double sigmaA, double sigmaB, double rc);
double energySystem(std::vector<std::vector<double>> positionArray, std::vector<double> radiusArray, double rc, double lengthCube);
double squareDistancePair(std::vector<double> positionA, std::vector<double> positionB, double lengthCube);
double energyParticle(int wu, int indexParticle, std::vector<double> positionParticle, std::vector<std::vector<double>> positionArray, std::vector<int> neighborIList, std::vector<double> radiusArray, double rc, double lengthCube);
std::vector<std::vector<int>> createNeighborList(std::vector<std::vector<double>> positionArray, double skin, double lengthCube);

#endif /* ENERGY_H_ */
