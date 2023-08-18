/*
 * random.h
 *
 *  Created on: 3 oct. 2022
 *      Author: Romain Simon
 */

#ifndef RANDOM_H_
#define RANDOM_H_

double randomDoubleGenerator(const double& min, const double& max);

int randomIntGenerator(const int& min, const int& max);

std::vector<double> randomVectorDoubleGenerator(const int& vectorSize, const double& min, const double& max);

std::vector<double> randomUnitSphereVectorGenerator(const int& vectorSize);

std::vector<double> randomRotationMatrixGenerator(const int& vectorSize);
#endif /* RANDOM_H_ */
