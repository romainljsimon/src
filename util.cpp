/*
 * util.cpp
 *
 *  Created on: 13 oct. 2022
 *      Author: rsimon
 */

#include <vector>

//Boundary conditions

std::vector<double> periodicBC(std::vector<double> positionParticle, double lengthCube)
/*
 *This function is an implementation of the periodic Boundary conditions.
 *If a particle gets out of the simulation box from one of the sides, it gets back in the box from the opposite side.
 */
{
	int positionParticleSize { static_cast<int>(positionParticle.size()) };
	for (int i = 0; i < positionParticleSize; i++)
	{
		if (positionParticle[i]< 0)
			positionParticle[i] += lengthCube;

		else if (positionParticle[i] > lengthCube)
			positionParticle[i] -= lengthCube;
	}
	return positionParticle;
}

// Vector Operators

double meanVector(std::vector<double> vec)
{
	double sum { 0 };
	for (auto const &elt: vec)
	{
		sum += elt;
	}
	return sum / static_cast<double>(vec.size());
}

std::vector<double> divideVectorByScalar(std::vector<double> vec, double scalar)
{
	for (int i = 0; i < static_cast<int>(vec.size()); i++)
	{
		vec[i] = vec[i] / scalar;
	}
	return vec;
}

std::vector<double> multiplyVectorByScalar(std::vector<double> vec, double scalar)
{
	for (int i = 0; i < static_cast<int>(vec.size()); i++)
	{
		vec[i] = vec[i] * scalar;
	}
	return vec;
}


std::vector<double> vectorNormalization(std::vector<double> vec)
{
	double mean { meanVector(vec) };
	return divideVectorByScalar(vec, mean);
}

std::vector<double> vectorSum(std::vector<double> vec1, std::vector<double> vec2)
{
	int vecSize{ static_cast<int>(vec1.size()) };
	for (int i = 0; i < vecSize; i++)
	{
		vec1[i] += vec2[i];
	}
	return vec1;
}

double getMaxVector(std::vector<double> vec)
{
	double max = vec[0];
	for (auto const &elt: vec)
	{
		max = (max < elt) ? elt: max;
	}
	return max;
}
//double meanMatrix(std::vector<std::vector<double>> mat)
//{
//	double sum{0};
//	for (auto const &vec: mat)
//	{
//		for (auto const &elt: vec)
//		{
//			sum += elt;
//		}
//
//	}
//	return sum / static_cast<double>(mat.size());
//}

std::vector<std::vector<double>> multiplyMatrixByScalar(std::vector<std::vector<double>> mat, double scalar)
{
	for (int i = 0; i < static_cast<int>(mat.size()); i++)
		{
			mat[i] = multiplyVectorByScalar(mat[i], scalar);
		}
	return mat;
}

double getMaxMatrix(std::vector<std::vector<double>> mat)
{
	double max { mat[0][0] };

	for (auto const &vec: mat)
	{
		double newMax = getMaxVector(vec);
		max = (newMax > max) ? newMax : max;
	}
	return max;
}

std::vector<std::vector<double>> rescaleMatrix(std::vector<std::vector<double>> mat, double rescaler)
{
	double max { getMaxMatrix(mat) };
	double ratio {rescaler / max};
	return multiplyMatrixByScalar(mat, ratio);
}
