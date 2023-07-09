/*
 * energy.h
 *
 *  Created on: 3 oct. 2022
 *      Author: Romain Simon
 */

#ifndef ENERGY_H_
#define ENERGY_H_


#include "NEIGHBORS/Neighbors.h"
#include "POTENTIALS/BondPotentials.h"
#include "POTENTIALS/PairPotentials.h"

double energyParticle(const int& indexParticle, const std::vector<double>& positionParticle,
                      const Particles& systemParticles, const std::vector<int>& neighborIList,
                      const PairPotentials& pairPotentials);

double energyParticle(const int& indexParticle, const std::vector<double>& positionParticle,
                      const Particles& systemParticles, const std::vector<int>& neighborIList,
                      const PairPotentials& pairPotentials, const int& indexSkip);

double energySystem(const Particles& systemParticles, const Neighbors& systemNeighbors,
                    const PairPotentials& pairPotentials);

double energyParticlePolymer (const int& indexParticle, const std::vector<double>& positionParticle,
                              const Particles& systemParticles, const std::vector<int>& neighborIList,
                              const PairPotentials& pairPotentials, const BondPotentials& bondPotentials,
                              const int& indexSkip);

double energyParticlePolymer (const int& indexParticle, const std::vector<double>& positionParticle,
                              const Particles& systemParticles, const std::vector<int>& neighborIList,
                              const PairPotentials& pairPotentials, const BondPotentials& bondPotentials);

double energySystemPolymer(const Particles& systemParticles, const Neighbors& systemNeighbors,
                           const PairPotentials& pairPotentials, const BondPotentials& bondPotentials);

#endif /* ENERGY_H_ */
