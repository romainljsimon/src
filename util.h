/*
 * util.h
 *
 *  Created on: 13 oct. 2022
 *      Author: rsimon
 */

#ifndef UTIL_H_
#define UTIL_H_

std::vector<double> periodicBC(std::vector<double> positionParticle, const double& lengthCube);

std::vector<double> divideVectorByScalar(std::vector<double> vec, const double& scalar);

std::vector<double> multiplyVectorByScalar(std::vector<double> vec, const double& scalar);

std::vector<double> vectorNormalization(std::vector<double> vec);

std::vector<double> vectorSum(std::vector<double> vec1, const std::vector<double>& vec2);

std::vector<double> vectorDiff(std::vector<double> vec1, std::vector<double> vec2);

double getMaxVector(const std::vector<double>& vec);

double getSquareNormVector(const std::vector<double>& vec);

std::vector<std::vector<double>> matrixSum(std::vector<std::vector<double>> mat1, const std::vector<std::vector<double>>& mat2);

std::vector<std::vector<double>> matrixSumWithVector(std::vector<std::vector<double>> mat1, const std::vector<double>& vec1);

std::vector<std::vector<double>> matrixDiffWithVector(std::vector<std::vector<double>> mat1, std::vector<double> vec1);

std::vector<std::vector<double>> multiplyMatrixByScalar(std::vector<std::vector<double>> mat, const double& scalar);

std::vector<std::vector<double>> rescaleMatrix(std::vector<std::vector<double>> mat, const double& rescaler);

std::vector<double> getSquareNormRowMatrix(std::vector<std::vector<double>> mat);

std::vector<double> meanColumnsMatrix(std::vector<std::vector<double>> mat);

std::vector<int> createSaveTime(const int& max, const int& linear_scalar, const float& log_scalar);

#endif /* UTIL_H_ */
