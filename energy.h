/*
 * energy.h
 *
 *  Created on: 3 oct. 2022
 *      Author: Romain Simon
 */

#ifndef ENERGY_H_
#define ENERGY_H_


double energySystem(const std::vector<std::vector<double>>& positionArray, const std::vector<double>& diameterArray,
                    const std::vector<std::vector<int>>& neighborList, const double& squareRc,
                    const double& lengthCube, const double& halfLengthCube);

double energyParticle(const int& indexParticle, const std::vector<double>& positionParticle,
		              const std::vector<std::vector<double>>& positionArray, const std::vector<int>& neighborIList,
					  const std::vector<double>& diameterArray, const double& squareRc, const double& lengthCube,
                      const double& halfLengthCube, const int& indexSkip = -1);

double energyParticlePolymer (const int& indexParticle, const std::vector<double>& positionParticle,
                              const std::vector<std::vector<double>>& positionArray,
                              const std::vector<int>& neighborIList, const std::vector<double>& diameterArray,
                              const std::vector<int>& bondsI, const double& squareRc, const double& lengthCube,
                              const double& halfLengthCube, const double& squareR0, const double& feneK,
                              const std::string& bondType, const int& indexSkip = -1);

double energySystemPolymer(const std::vector<std::vector<double>>& positionArray,
                           const std::vector<double>& diameterArray,
						   const std::vector<std::vector<int>>& bondsMatrix,
                           const std::vector<std::vector<int>>& neighborList,
						   const double& squareRc, const double& lengthCube, const double& halfLengthCube,
                           const double& squareR0, const double& feneK, const std::string& bond_Type);

#endif /* ENERGY_H_ */
