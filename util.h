/*
 * util.h
 *
 *  Created on: 13 oct. 2022
 *      Author: Romain Simon
 */

#ifndef UTIL_H_
#define UTIL_H_

#include <vector>
#include <numeric>
#include <algorithm>
#include <math.h>

/***
double squareDistancePair(const std::vector<double>& positionA,  const std::vector<double>& positionB,
                          const double& lengthCube, const double& halfLengthCube);
***/
std::vector<double> periodicBC(std::vector<double> positionParticle, const double& lengthCube);

std::vector<double> divideVectorByScalar(std::vector<double> vec, const double& scalar);

std::vector<double> multiplyVectorByScalar(std::vector<double> vec, const double& scalar);


std::vector<double> vectorSum(const std::vector<double>& vec1, const std::vector<double>& vec2);

std::vector<double> vectorDiff(std::vector<double> vec1, std::vector<double> vec2);

double getMaxVector(const std::vector<double>& vec);

double getSquareNormVector(const std::vector<double>& vec);

std::vector<std::vector<double>> matrixSum(std::vector<std::vector<double>> mat1, const std::vector<std::vector<double>>& mat2);

std::vector<std::vector<double>> matrixSumWithVector(std::vector<std::vector<double>> mat1, const std::vector<double>& vec1);

std::vector<std::vector<double>> matrixDiffWithVector(std::vector<std::vector<double>> mat1, std::vector<double> vec1);

std::vector<std::vector<double>> multiplyMatrixByScalar(std::vector<std::vector<double>> mat, const double& scalar);

std::vector<std::vector<double>> rescaleMatrix(std::vector<std::vector<double>> mat, const double& rescale);

template<typename InputIt>
double getSquareNormVector(const InputIt& vecItBegin, const InputIt& vecItEnd)
{
    return std::inner_product(vecItBegin, vecItEnd, vecItBegin, 0.);
}

template<typename MatInputIt, typename VecInputIt>
std::vector<double> matrixMultiplication(MatInputIt matItBegin, int lenMat, VecInputIt vecItBegin, int lenVec)
{
    const int nRows {lenMat / lenVec};
    std::vector<double> out;
    for (int i=0; i < nRows; ++i)
    {
        std::vector<double> rowMultiplication(lenVec);
        std::transform(vecItBegin, vecItBegin + lenVec,
                       matItBegin + i * lenVec, rowMultiplication.begin(), std::multiplies<>());
        out.push_back(std::accumulate(rowMultiplication.begin(), rowMultiplication.begin() + lenVec, 0.));
    }
    return out;
}

template<typename InputIt>
std::vector<double> getSquareNormRowMatrix(InputIt vecItBegin, const int& lenColumn, const int& lenRow)
{
    std::vector<double> squareNormVector;
    auto vecItEnd {vecItBegin + lenRow};
    for (int i = 0; i < lenColumn; i++)
    {
        squareNormVector.push_back(getSquareNormVector( vecItBegin, vecItEnd));
        vecItBegin = vecItEnd;
        vecItEnd = vecItBegin + lenRow;
    }
    return squareNormVector;
}
template<typename InputIt>
void vectorNormalization(InputIt vecItBegin, InputIt vecItEnd)
{
    double squareNorm {getSquareNormVector(vecItBegin, vecItEnd)};
    double norm {pow(squareNorm, 1./2)};
    std::for_each(vecItBegin, vecItEnd, [&norm](double &n) { n / norm; });
}

std::vector<double> meanColumnsMatrix(std::vector<std::vector<double>> mat);

std::vector<int> createSaveTime(const int& max, const int& linear_scalar, const float& log_scalar);

#endif /* UTIL_H_ */
