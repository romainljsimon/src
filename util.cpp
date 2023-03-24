/*
 * util.cpp
 *
 *  Created on: 13 oct. 2022
 *      Author: rsimon
 */

#include <vector>
#include <functional>
#include <numeric>

//Boundary conditions

std::vector<double> periodicBC(std::vector<double> positionParticle, const double& lengthCube)
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

std::vector<double> divideVectorByScalar(std::vector<double> vec, const double& scalar)
{
	for (int i = 0; i < static_cast<int>(vec.size()); i++)
	{
		vec[i] = vec[i] / scalar;
	}
	return vec;
}

std::vector<double> multiplyVectorByScalar(std::vector<double> vec, const double& scalar)
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

std::vector<double> vectorSum(std::vector<double> vec1, const std::vector<double>& vec2)
{
	std::transform (vec1.begin(), vec1.end(), vec2.begin(), vec1.begin(), std::plus<double>());
	return vec1;
}

std::vector<double> vectorDiff(std::vector<double> vec1, std::vector<double> vec2)
{
	std::vector<double> minus_vec2 = multiplyVectorByScalar(vec2, -1.);
	return vectorSum(vec1, minus_vec2);
}

double getMaxVector(const std::vector<double>& vec)
{
	double max = vec[0];
	for (auto const &elt: vec)
	{
		max = (max < elt) ? elt: max;
	}
	return max;
}

double getSquareNormVector(const std::vector<double>& vec)
{
	return std::inner_product(vec.begin(), vec.end(), vec.begin(), 0.);
}

// matrix operations


std::vector<std::vector<double>> matrixSum(std::vector<std::vector<double>> mat1, const std::vector<std::vector<double>>& mat2)
{
	int matSize{ static_cast<int>(mat1.size()) };

	for (int i = 0; i < matSize; i++)
	{
			mat1[i] = vectorSum(mat1[i], mat2[i]);
	}
	return mat1;

}

std::vector<std::vector<double>> matrixSumWithVector(std::vector<std::vector<double>> mat1, const std::vector<double>& vec1)
{
	int matSize{ static_cast<int>(mat1.size()) };

	for (int i = 0; i < matSize; i++)
	{
			mat1[i] = vectorSum(mat1[i], vec1);
	}
	return mat1;

}

std::vector<std::vector<double>> matrixDiffWithVector(std::vector<std::vector<double>> mat1, std::vector<double> vec1)
{
	int matSize{ static_cast<int>(mat1.size()) };

	for (int i = 0; i < matSize; i++)
	{
			mat1[i] = vectorDiff(mat1[i], vec1);
	}
	return mat1;
}

std::vector<std::vector<double>> multiplyMatrixByScalar(std::vector<std::vector<double>> mat, const double& scalar)
{
	for (int i = 0; i < static_cast<int>(mat.size()); i++)
		{
			mat[i] = multiplyVectorByScalar(mat[i], scalar);
		}
	return mat;
}

double getMaxMatrix(const std::vector<std::vector<double>>& mat)
{
	double max { mat[0][0] };

	for (auto const &vec: mat)
	{
		double newMax = getMaxVector(vec);
		max = (newMax > max) ? newMax : max;
	}
	return max;
}

std::vector<std::vector<double>> rescaleMatrix(std::vector<std::vector<double>> mat, const double& rescaler)
{
	double max { getMaxMatrix(mat) };
	double ratio {rescaler / max};
	return multiplyMatrixByScalar(mat, ratio);
}

std::vector<double> getSquareNormRowMatrix(std::vector<std::vector<double>> mat)
{
	int matSize { static_cast<int>(mat.size()) };
	std::vector<double> squareNormVector ( matSize );
	for (int i = 0; i < matSize; i++)
	{
		squareNormVector[i] = getSquareNormVector( mat[i] );
	}
	return squareNormVector;
}

std::vector<double> meanColumnsMatrix(std::vector<std::vector<double>> mat)
{
	int matSize { static_cast<int>(mat.size()) };
	std::vector<double> mean(3, 0);

	for (int i = 0; i < matSize; i++)
	{
		mean = vectorSum(mean, mat[i]);
	}
	return divideVectorByScalar(mean, matSize);
}

std::vector<int> createSaveTime(const int& max, const int& linear_scalar, const float& log_scalar)
{
	std::vector<int> timeStepArray;

	for (int j = 0; j < max; j += linear_scalar)
	{
		timeStepArray.push_back ( j ) ;
		timeStepArray.push_back ( j + 1 );


		for (int i = static_cast<int>(log_scalar) +1; i < linear_scalar; i = static_cast<int>(i * log_scalar) +1)
		{
			timeStepArray.push_back ( j + i );

		}
	}
	timeStepArray.push_back ( max ) ;
	return timeStepArray;

}
