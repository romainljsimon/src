/*
 * random.h
 *
 *  Created on: 3 oct. 2022
 *      Author: romainsimon
 */

#ifndef RANDOM_H_
#define RANDOM_H_

double randomDoubleGenerator(const double& min, const double& max);

int randomIntGenerator(const int& min, const int& max);

std::vector<double> randomVectorDoubleGenerator(const int& vectorSize, const double& min, const double& max);
#endif /* RANDOM_H_ */
