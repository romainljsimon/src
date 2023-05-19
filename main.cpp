//============================================================================
// Name        : main.cpp
// Author      : Romain Simon
// Version     :
// Copyright   : Your copyright notice
// Description : main function of the swapMC project in C++, Ansi-style
//============================================================================

#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <chrono>
#include "readSaveFile.h"
#include "util.h"
#include "MonteCarlo.h"
#include <input/Parameter.h>


int main()
{
	std::string folderPath ( "." );


	//Opening input variables file
    param::Parameter param(folderPath + "/inputVar.txt" );
	posRad initPosRad = readXYZ ( folderPath + "/initPosition.xyz");

    // double density =  param.get_double("density");
	// std::filesystem::create_directory (folderPath + "/outXYZ" );
	// std::filesystem::create_directory (folderPath + "/disp" );


    //Opening position file
	initPosRad.radVector = vectorNormalization(initPosRad.radVector);
    // We define the lengths cale of the system as <sigma>. radVector is now the diameterVector
	// initPosRad.radVector = divideVectorByScalar(initPosRad.radVector, 2);

	//assert((initVar.rc > 1.) && "The cut-off radius has to be superior to 1.");
	//assert((initVar.temp > 0.) && "The temperature has to be superior to 0.");
	//assert((initVar.density > 0.) && "The number density has to be superior to 0.");
	//assert((initVar.rBox > 0.) && "The length of the translation box has to be superior to 0.");
	//assert((lengthCube > 1.) && "The length of the system has to be bigger than a particle");

    //initPosRad.posMatrix = rescaleMatrix(initPosRad.posMatrix, lengthCube);
	//std::cout << squareDistancePairTest();
	//randomGeneratorTest();
    // mcTotal(initPosRad.posMatrix, initPosRad.radVector,  rc, lengthCube, temp, rBox, number);
	//mcTotal(initPosRad.posMatrix, initPosRad.radVector,  rc, lengthCube, temp, rBox, number);
	std::clock_t c_start = std::clock();
	auto t_start = std::chrono::high_resolution_clock::now();
    std::string key { "simType" };
    std::string simType { param.get_string(key) };

    if (simType == "polymer")
    {
        std::vector<std::vector<int>> bondsMatrix = readBondsTXT(folderPath + "/bonds.txt");

        MonteCarlo system {param, initPosRad.posMatrix, initPosRad.radVector,
                          initPosRad.moleculeType, bondsMatrix, folderPath};

        system.mcTotal();
    }
    else if (simType == "atomic")
    {
        MonteCarlo system{param,initPosRad.posMatrix,initPosRad.radVector,
                          initPosRad.moleculeType, folderPath};

        system.mcTotal();
    }
    else
    {
        std::cout << "The simType variable is not 'atomic' or 'polymer' \n";
        std::cout << "Please define the simulation type as one of the two\n";
    }

	std::clock_t c_end = std::clock();
	auto t_end = std::chrono::high_resolution_clock::now();
	auto cpuTime {(static_cast<double>(c_end - c_start)) / CLOCKS_PER_SEC};
	auto wallTime {std::chrono::duration<double, std::milli>(t_end - t_start).count() / 1000.};

	std::cout << "CPU time used: " << cpuTime << " seconds; "
			  << cpuTime / 60. << " minutes; " << cpuTime / 3600. << " hours\n";
	std::cout << "Wall clock time passed: " << wallTime << " seconds; "
			  << wallTime / 60. << " minutes; " << wallTime / 3600. << " hours\n";

	return 0;
}
