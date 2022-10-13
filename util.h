/*
 * util.h
 *
 *  Created on: 13 oct. 2022
 *      Author: rsimon
 */

#ifndef UTIL_H_
#define UTIL_H_

std::vector<double> divideVectorByScalar(std::vector<double> vec, double scalar);
std::vector<double> vectorNormalization(std::vector<double> vec);
std::vector<std::vector<double>> multiplyMatrixByScalar(std::vector<std::vector<double>> mat, double scalar);
std::vector<std::vector<double>> rescaleMatrix(std::vector<std::vector<double>> mat, double rescaler);

#endif /* UTIL_H_ */
