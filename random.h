/*
 * random.h
 *
 *  Created on: 3 oct. 2022
 *      Author: romainsimon
 */

#ifndef RANDOM_H_
#define RANDOM_H_

double randomDoubleGenerator(double min, double max);
int randomIntGenerator(int min, int max);
std::vector<double> randomVectorDoubleGenerator(int vectorSize, double min, double max);

#endif /* RANDOM_H_ */
