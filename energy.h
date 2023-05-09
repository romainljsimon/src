/*
 * energy.h
 *
 *  Created on: 3 oct. 2022
 *      Author: romainsimon
 */

#ifndef ENERGY_H_
#define ENERGY_H_

double ljPotential(const double& squareDistance, const double& sigmaA, const double& sigmaB, const double& squareRc,
                   const double& shift);

double squareDistancePair(const std::vector<double>& positionA,  const std::vector<double>& positionB, const double& lengthCube);

double energySystem(const std::vector<std::vector<double>>& positionArray, const std::vector<double>& diameterArray,
                    const std::vector<std::vector<int>>& neighborList, const double& squareRc,
                    const double& lengthCube);

double energyParticle(const int& indexParticle, const std::vector<double>& positionParticle,
                      const std::vector<std::vector<double>>& positionArray, const std::vector<int>& neighborIList,
                      const std::vector<double>& diameterArray, const double& squareRc, const double& lengthCube,
                      const int& indexSkip = -1);

double energyParticlePolymer (const int& indexParticle, const std::vector<double>& positionParticle,
                              const std::vector<std::vector<double>>& positionArray,
                              const std::vector<int>& neighborIList, const std::vector<double>& diameterArray,
                              const std::vector<int>& bondsI, const double& squareRc, const double& lengthCube,
                              const double& squareR0, const double& feneK, const int& indexSkip = -1);

double energySystemPolymer(const std::vector<std::vector<double>>& positionArray,
                           const std::vector<double>& diameterArray,
						   const std::vector<std::vector<int>>& bondsMatrix,
                           const std::vector<std::vector<int>>& neighborList,
						   const double& squareRc, const double& lengthCube,
                           const double& squareR0, const double& feneK);

double pressureSystem(double temp, std::vector<std::vector<double>> positionArray, std::vector<double> diameterArray,
                      double squareRc, double lengthCube);

double pressureParticle(double temp, int indexParticle, std::vector<double> positionParticle,
                        std::vector<std::vector<double>> positionArray, std::vector<int> neighborIList,
                        std::vector<double> diameterArray, double squareRc, double lengthCube);




#endif /* ENERGY_H_ */
