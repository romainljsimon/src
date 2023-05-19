/*
 * pressure.h
 *
 *  Created on: 11 apr. 2023
 *      Author: Romain Simon
 */

#ifndef PRESSURE_H
#define PRESSURE_H


double pressureParticle(const double& temp, int indexParticle, const std::vector<double>& positionParticle,
                        const std::vector<std::vector<double>>& positionArray, const std::vector<int>& neighborIList,
                        const std::vector<double>& diameterArray, const double& squareRc, const double& lengthCube);

double pressureSystem(const double& temp, const std::vector<std::vector<double>>& positionArray,
                      const std::vector<double>& diameterArray, const double& squareRc, const double& lengthCube);

#endif /* PRESSURE_H */
