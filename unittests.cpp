/*
 * unittests.cpp
 *
 *  Created on: 6 oct. 2022
 *      Author: romainsimon
 */

#include <iostream>
#include <vector>
#include "energy.h"
#include "random.h"


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

	for (int i = 0; i < static_cast<int>(distanceArray.size()); i++)
	{
		double squareDistance {squareDistancePair(positionA, positionArrayB[i], lengthCube)};

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
	constexpr int nrolls     = 10000;  // number of experiments
	constexpr int nstars     = 95;     // maximum number of stars to distribute
	constexpr int nintervals = 6;      // number of intervals


	int p[nintervals]={};

	for (int i=0; i < nrolls; ++i)
	{
		int number = randomIntGenerator(0, 5);
	    ++p[number];
	}

	std::cout << "uniform_discrete_distribution (0, 6):" << '\n';
	std::cout << std::fixed; std::cout.precision(1);

	for (int i=0; i < nintervals; ++i)
	{
		std::cout << i << ": ";
	    std::cout << std::string(p[i]*nstars/nrolls,'*') << std::endl;
	}

}

void randomDoubleGeneratorTest()
{
  constexpr int nrolls=10000;  // number of experiments
  constexpr int nstars=95;     // maximum number of stars to distribute
  constexpr int nintervals=10; // number of intervals

  int p[nintervals]={};

  for (int i=0; i<nrolls; ++i) {
    double number = randomDoubleGenerator(0., 1.);
    ++p[int(nintervals*number)];
  }

  std::cout << "uniform_real_distribution (0.0,1.0):" << '\n';
  std::cout << std::fixed; std::cout.precision(1);

  for (int i=0; i < nintervals; ++i) {
    std::cout << float(i) / nintervals << "-" << float(i+1) / nintervals << ": ";
    std::cout << std::string(p[i] * nstars / nrolls,'*') << '\n';
  }

}

void randomGeneratorTest()
{
	randomDoubleGeneratorTest();
	randomIntGeneratorTest();
}



