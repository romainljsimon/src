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
double energyParticle(int indexParticle, std::vector<double> positionParticle, std::vector<std::vector<double>> positionArray, std::vector<double> radiusArray, double rc, double lengthCube);

#endif /* ENERGY_H_ */
