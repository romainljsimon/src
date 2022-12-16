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
std::vector<double> vectorSum(std::vector<double> vec1, std::vector<double> vec2);
double getMaxVector(std::vector<double> vec);
std::vector<std::vector<double>> matrixSum(std::vector<std::vector<double>> mat1, std::vector<std::vector<double>> mat2);
std::vector<std::vector<double>> multiplyMatrixByScalar(std::vector<std::vector<double>> mat, double scalar);
std::vector<std::vector<double>> rescaleMatrix(std::vector<std::vector<double>> mat, double rescaler);
std::vector<double> vectorSum(std::vector<double> vec1, std::vector<double> vec2);
std::vector<double> periodicBC(std::vector<double> positionParticle, double lengthCube);
std::vector<double> getSquareNormRowMatrix(std::vector<std::vector<double>> mat);
std::vector<double> meanColumnsMatrix(std::vector<std::vector<double>> mat);
std::vector<double> vectorDiff(std::vector<double> vec1, std::vector<double> vec2);
std::vector<std::vector<double>> matrixSumWithVector(std::vector<std::vector<double>> mat1, std::vector<double> vec1);
std::vector<std::vector<double>> matrixDiffWithVector(std::vector<std::vector<double>> mat1, std::vector<double> vec1);
double getSquareNormVector(std::vector<double> vec);

#endif /* UTIL_H_ */
