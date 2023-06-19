/*
 * unittests.cpp
 *
 *  Created on: 6 oct. 2022
 *      Author: Romain Simon
 */

#include <iostream>
#include <vector>
#include "random.h"
#include "util.h"

int squareDistancePairTest()
{
	std::vector<std::vector<double>> positionArrayB { { 0., 0., 0. }, { 1., 1., 1. }, { 1., 0., 0. }, { 0., 1., 0. },
														{ 0., 0., 1. }, { 1., 0., 1. }, { 0., 1., 1. }, { 1., 1., 0. },
                                                        { 0.5, 0., 0. }, { 0., 0.5, 0. }, { 0., 0., 0.5 }, { 0., 0.5, 0.5 },
                                                            { 0.5, 0., 0.5 }, { 0.5, 0.5, 0. }, { 0.5, 0.5, 0.5 }};

	std::vector<double> positionA { 0., 0., 0. };
	std::vector<double> distanceArray { 0., 0., 0., 0.,
										0., 0., 0., 0.,
	 	 	 	 	 	 	 	 	 	0.25, 0.25, 0.25, 0.5,
										0.5, 0.5, 0.75};

	constexpr double lengthCube { 1. };
    constexpr double halfLengthCube { 1. };

	for (int i = 0; i < static_cast<int>(distanceArray.size()); i++)
	{
		double squareDistance {squareDistancePair(positionA, positionArrayB[i], lengthCube, halfLengthCube)};

		if (squareDistance != distanceArray[i])
		{
			std::cout << i << " " << squareDistance << " " << distanceArray[i] << " "<< "\n";
			return 1;
		}

	}
	return 0;
}

void randomIntGeneratorTest()
{
	constexpr int nRolls     = 10000;  // number of experiments
	constexpr int nStars     = 95;     // maximum number of stars to distribute
	constexpr int nIntervals = 6;      // number of intervals


	int p[nIntervals]={};

	for (int i=0; i < nRolls; ++i)
	{
		int number = randomIntGenerator(0, 5);
	    ++p[number];
	}

	std::cout << "uniform_discrete_distribution (0, 6):" << '\n';
	std::cout << std::fixed; std::cout.precision(1);

	for (int i=0; i < nIntervals; ++i)
	{
		std::cout << i << ": ";
	    std::cout << std::string(p[i]*nStars/nRolls,'*') << std::endl;
	}

}

void randomDoubleGeneratorTest()
{
  constexpr int nRolls=10000;  // number of experiments
  constexpr int nStars=95;     // maximum number of stars to distribute
  constexpr int nIntervals=10; // number of intervals

  int p[nIntervals]={};

  for (int i=0; i<nRolls; ++i) {
    double number = randomDoubleGenerator(0., 1.);
    ++p[int(nIntervals*number)];
  }

  std::cout << "uniform_real_distribution (0.0,1.0):" << '\n';
  std::cout << std::fixed; std::cout.precision(1);

  for (int i=0; i < nIntervals; ++i) {
    std::cout << float(i) / nIntervals << "-" << float(i+1) / nIntervals << ": ";
    std::cout << std::string(p[i] * nStars / nRolls,'*') << '\n';
  }

}

void randomGeneratorTest()
{
	randomDoubleGeneratorTest();
	randomIntGeneratorTest();
}



