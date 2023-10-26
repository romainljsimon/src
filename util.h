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

/***
double squareDistancePair(const std::vector<double>& positionA,  const std::vector<double>& positionB,
                          const double& lengthCube, const double& halfLengthCube);
***/
double innerProduct(const std::vector<double>& vec1, const std::vector<double>& vec2);

double cosAngleVectors(const std::vector<double>& vec1, const std::vector<double>& vec2);

std::vector<double> divideVectorByScalar(std::vector<double> vec, const double& scalar);

std::vector<double> multiplyVectorByScalar(std::vector<double> vec, const double& scalar);

std::vector<double> vectorNormalization(const std::vector<double>& vec);

std::vector<double> vectorSum(const std::vector<double>& vec1, const std::vector<double>& vec2);

std::vector<double> vectorDiff(const std::vector<double>& vec1, const std::vector<double>& vec2, int nDims);

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
void periodicVector(InputIt vecItBegin, const int& lenVec, const double& lengthCube)
{
    std::for_each(vecItBegin, vecItBegin + lenVec, [&lengthCube](auto &n) { n = (n < -lengthCube/2) ? n + lengthCube: n;});
    std::for_each(vecItBegin, vecItBegin + lenVec, [&lengthCube](auto &n) { n = (n > lengthCube/2) ? n - lengthCube: n;});
}

template<typename InputIt>
void sumVector(InputIt vecItBegin1, InputIt vecItBegin2, const int& lenVec)
{
    std::transform (vecItBegin1, vecItBegin1+lenVec, vecItBegin2, vecItBegin1, std::plus<>());
}
template<typename InputIt>
void divVector(InputIt vecItBegin, const int& lenVec, const double& fac)
{
    std::for_each(vecItBegin, vecItBegin + lenVec, [&fac](auto &n) { n /= fac;});
}

template<typename InputIt>
std::vector<double> crossProduct(InputIt vecItBegin1, InputIt vecItBegin2, const int& lenVec)
{
    std::vector<double> cross (lenVec);
    cross[0] = (*(vecItBegin1 + 1)) * (*(vecItBegin2 + 2)) - (*(vecItBegin1 + 2)) * (*(vecItBegin2 + 1));
    cross[1] = (*(vecItBegin1 + 2)) * (*(vecItBegin2)) - (*(vecItBegin1)) * (*(vecItBegin2 + 2));
    cross[2] = (*(vecItBegin1)) * (*(vecItBegin2 + 1)) - (*(vecItBegin1 + 1)) * (*(vecItBegin2));
    return cross;
}

template<typename InputIt>
double traceMatrix(InputIt vecItBegin, const int& lenRowCol)
{
    double trace {0.};
    int lenMat { lenRowCol * lenRowCol};
    for (auto it=vecItBegin; it < (vecItBegin + lenMat); it+=(lenRowCol+1))
    {
        trace += *it;
    }
    return trace;
}


std::vector<double> meanColumnsMatrix(std::vector<std::vector<double>> mat);

std::vector<int> createSaveTime(const int& max, const int& linear_scalar, const float& log_scalar);

#endif /* UTIL_H_ */
