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
	posRad initPosRad = readXYZ("/home/rsimon/eclipse-workspace/swapMC/input/inputt4/randomt4v1.xyz");
    //Opening the file
	initPosRad.radVector = vectorNormalization(initPosRad.radVector);     // We define the lengthscale of the system as <sigma>
	initPosRad.radVector = divideVectorByScalar(initPosRad.radVector, 2);
	double rc { 2.5 };
	double temp { 4. };
	double density { 1. };
	double rbox{ 0.3 };
	double lengthCube {pow (static_cast<double>(initPosRad.radVector.size()) / density, 1./3) };
	initPosRad.posMatrix = rescaleMatrix(initPosRad.posMatrix, lengthCube);
	//std::cout << squareDistancePairTest();
	//randomGeneratorTest();
    mcTotal(initPosRad.posMatrix, initPosRad.radVector,  rc, lengthCube, temp, rbox, number);
	return 0;
}
