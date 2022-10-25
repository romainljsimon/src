//============================================================================
// Name        : main.cpp
// Author      : Romain Simon
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include <vector>
#include <cassert>
#include "unittests.h"
#include "readSaveFile.h"
#include "random.h"
#include "energy.h"
#include "mc.h"
#include "util.h"


int main()
{
	std::string number;
	std::cout << "Number of the simulation: ";
	std::cin >> number;
	posRad initPosRad = readXYZ("/Users/romainsimon/eclipse-workspace/swapMC/input/inputt4/randomv1t4.xyz");
    //Opening the file
	initPosRad.radVector = vectorNormalization(initPosRad.radVector);     // We define the lengthscale of the system as <sigma>
	initPosRad.radVector = divideVectorByScalar(initPosRad.radVector, 2);
	double rc { 2.5 };
	double temp { 4. };
	double density { 1. };
	double rbox{ 0.3 };
	double lengthCube {pow (static_cast<double>(initPosRad.radVector.size()) / density, 1./3) };
	assert((rc > 1.) && "The cut-off radius has to be superior to 1.");
	assert((temp > O.) && "The temperature has to be superior to 0.");
	assert((density > 0.) && "The number density has to be superior to 0.");
	assert((rbox > 0.) && "The length of the translation box has to be superior to 0.");
	assert((lengthCube > 1.) && "The length of the system has to be bigger than a particle");

	initPosRad.posMatrix = rescaleMatrix(initPosRad.posMatrix, lengthCube);
	//std::cout << squareDistancePairTest();
	//randomGeneratorTest();
    mcTotal(initPosRad.posMatrix, initPosRad.radVector,  rc, lengthCube, temp, rbox, number);
	return 0;
}
