/*
 * pressure.cpp
 *
 *  Created on: 11 apr. 2023
 *      Author: Romain Simon
 */

#include <vector>
#include <cmath>
#include "util.h"
/***
double rForce(double squareDistance, double sigmaA, double sigmaB, double squareRc)
{
    double squareSigma { std::pow (((sigmaA + sigmaB) / 2. ), 2.) };
    if (squareDistance > squareRc * squareSigma) // We consider a cut-off radius which is the threshold maximum distance of interaction between two particles
    {
        return 0.;
    }
    else
    {
        return 24. * (std::pow (( squareSigma / squareDistance), 3. ) - 2. * std::pow (( squareSigma / squareDistance ), 6.));
    }
}


double pressureParticle(const double& temp, const int indexParticle, const std::vector<double>& positionParticle,
                        const std::vector<std::vector<double>>& positionArray, const std::vector<int>& neighborIList,
                        const std::vector<double>& diameterArray, const double& squareRc, const double& lengthCube,
                        const double& halfLengthCube)
{
    double pressure { 0. };
    double particleDiameter = diameterArray[indexParticle];
    int neighborIListSize { static_cast<int>(neighborIList.size()) };
    int nParticles {static_cast<int>(positionArray.size())};

    for (int i = 0; i < neighborIListSize; i++)
    {
        int realIndex = neighborIList[i];
        double squareDistance { squareDistancePair (positionParticle, positionArray[realIndex], lengthCube, halfLengthCube) };

        if (realIndex == indexParticle)
        {
            continue;
        }

        pressure += rForce(squareDistance, particleDiameter, diameterArray[realIndex], squareRc);
    }
    return  - pressure / (3. * temp *  nParticles);
}


double pressureSystem(const double& temp, const std::vector<std::vector<double>>& positionArray,
                      const std::vector<double>& diameterArray, const double& squareRc, const double& lengthCube, const double& halfLengthCube)
{
    double sum { 0 };
    int positionArraySize {static_cast<int>(positionArray.size())};
    for (int i = 0; i < positionArraySize - 1; i++)
    {
        for (int j = i + 1; j < positionArraySize; j++) //inner loop for columns
        {
            double squareDistance { squareDistancePair (positionArray[i], positionArray[j], lengthCube, halfLengthCube)};
            sum += rForce(squareDistance, 2. * diameterArray[i], 2. * diameterArray[j], squareRc);
        }

    }
    return 1. - sum / (3. *  positionArraySize * temp);
}
***/
