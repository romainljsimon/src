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
#include "MonteCarlo.h"
#include <filesystem>
#include <chrono>

class Timer
{
private:
    // Type aliases to make accessing nested type easier
    using Clock = std::chrono::steady_clock;
    using Second = std::chrono::duration<double, std::ratio<1> >;

    std::chrono::time_point<Clock> m_beg{ Clock::now() };

public:

    void reset()
    {
        m_beg = Clock::now();
    }

    double elapsed() const
    {
        return std::chrono::duration_cast<Second>(Clock::now() - m_beg).count();
    }
};

int main()
{
	std::string folderPath ( std::filesystem::current_path() );


	//Opening input variables file

	inputVar initVar = readInput ( folderPath + "/inputVar.txt" );
	posRad initPosRad = readXYZ ( folderPath + "/initPosition.xyz", initVar.simulationMol);
	std::filesystem::create_directory (folderPath + "/outXYZ" );
	std::filesystem::create_directory (folderPath + "/disp" );


    //Opening position file
	initPosRad.radVector = vectorNormalization(initPosRad.radVector);     // We define the lengthscale of the system as <sigma>
	initPosRad.radVector = divideVectorByScalar(initPosRad.radVector, 2);

	//

	double lengthCube {pow (static_cast<double>(initPosRad.radVector.size()) / initVar.density, 1./3) };

	//assert((initVar.rc > 1.) && "The cut-off radius has to be superior to 1.");
	//assert((initVar.temp > 0.) && "The temperature has to be superior to 0.");
	//assert((initVar.density > 0.) && "The number density has to be superior to 0.");
	//assert((initVar.rbox > 0.) && "The length of the translation box has to be superior to 0.");
	//assert((lengthCube > 1.) && "The length of the system has to be bigger than a particle");

    // initPosRad.posMatrix = rescaleMatrix(initPosRad.posMatrix, lengthCube);
	//std::cout << squareDistancePairTest();
	//randomGeneratorTest();
    // mcTotal(initPosRad.posMatrix, initPosRad.radVector,  rc, lengthCube, temp, rbox, number);
	//mcTotal(initPosRad.posMatrix, initPosRad.radVector,  rc, lengthCube, temp, rbox, number);
	Timer t;
	MonteCarlo system { initVar.simulationMol, initPosRad.posMatrix, initPosRad.radVector,
						initPosRad.moleculeType,initVar.rc, lengthCube,
						initVar.temp, initVar.rbox, initVar.rskin, initVar.neighUpdate,
						folderPath, initVar.neighborMethod, initVar.timeSteps, initVar.r0, initVar.feneK};
	system.mcTotal();

    std::cout << "Time taken: " << t.elapsed() << " seconds\n";

	std::ofstream outTime;
	outTime.open(folderPath + "/time.txt", std::ios_base::app);
	outTime << t.elapsed();
	outTime << "\n";
	outTime.close();

	return 0;
}
