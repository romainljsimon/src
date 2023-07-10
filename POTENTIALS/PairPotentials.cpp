

#include <vector>
#include <cmath>
#include "PairPotentials.h"


int PairPotentials::getIndexIJ(const int& i, const int& j) const
{
    int indexI{ i };
    int indexJ{ j };

    if (i > j)
    {
        indexI = j;
        indexJ = i;
    }
    const int indexIJ {indexJ - indexI + m_nParticleTypes * ( indexI - 1) - ((indexI - 1) * (indexI-2)) / 2};
    return indexIJ * 4;
}

int PairPotentials::getParticleTypes() const
{
    return m_nParticleTypes;
}

double PairPotentials::getSquareRcIJ(const int& i, const int& j) const
{
    const int rcIndex {getIndexIJ(i, j) + 2};
    return m_pairPotentials[rcIndex];
}
/***
std::vector<double> PairPotentials::getPotentialsIJ(const int& i, const int& j, const int& k) const
{
    const int indexIJ { getIndexIJ(i, j) };
    return m_pairPotentials[indexIJ];
}
std::vector<double> PairPotentials::getPotentialsIJ(const int& i, const int& j) const
{
    const int indexIJ { getIndexIJ(i, j) };
    return m_pairPotentials[indexIJ];
}
 ***/
/*******************************************************************************
* This function calculates the Lennard-Jones potential energy between two
* particles seperated by a distance whose square is equal to squareDistance.
* We consider a cut off whose square is equal to squareRc times the mean diameter
* of the two particles considered. The potential is equal to 0 for distances
* greater that the cut off.
*
* @param squareDistance Square of the distance separating the two considered
*                       particles
*        sigmaA Particle A's diameter
*        sigmaB Particle B's diameter
*        squareRc Square of the cut off radius
*        shift Shifts the Lennard-Jones potential by a constant value
*
* @return energy
******************************************************************************/
double PairPotentials::ljPairEnergy(const double& squareDistance, const int& typeI, const int& typeJ) const
{
    //const std::vector<double>& potentialsIJ (getPotentialsIJ(typeI, typeJ));
    const int indexIJ { getIndexIJ(typeI, typeJ) };
    const double& rcSquareIJ {m_pairPotentials[indexIJ + 2]};

    if (squareDistance > rcSquareIJ)
    {
        return 0.;
    }

    else
    {
        const double& epsilonIJ { m_pairPotentials[indexIJ + 0]};
        const double& squareSigmaIJ { m_pairPotentials[indexIJ + 1]};
        const double& shiftIJ {m_pairPotentials[indexIJ + 3] };
        const double rapSquare { squareSigmaIJ / squareDistance };
        const double rapSix {rapSquare * rapSquare * rapSquare};
        return 4. * epsilonIJ * rapSix * ( rapSix - 1.) + 4. * epsilonIJ * shiftIJ;
    }
}

/*******************************************************************************
 * This function calculates the potential energy of one particle considering
 * that particles are Lennard-Jones particles.
 *
 * @param indexParticle Considered particle's index in the positionArray and
 *                      typeArray arrays.
 *        positionParticle Considered particle's position.
 * 	      positionArray array of particle positions.
 * 	      neighborIList Considered particle's neighbor list.
 *        typeArray array of particle diameters.
 *        squareRc square of the cut off radius.
 *        lengthCube simulation box length.
 * @return Particle's total energy.
 ******************************************************************************/
double PairPotentials::energyPairParticle(const int& indexParticle, const std::vector<double>& positionParticle,
                                          const Particles& systemParticles, const std::vector<int>& neighborIList,
                                          const int& indexSkip = -1) const
{
    double energy { 0. };
    const int particleType = systemParticles.getParticleTypeI(indexParticle);
    const int neighborIListSize { static_cast<int>(neighborIList.size()) };


    for (int i = 0; i < neighborIListSize; i++)
    {
        int realIndex = neighborIList[i];

        if (realIndex != indexSkip)
        {
            //if (realIndex == indexParticle)
            //{
            //    continue;
            //}

            const double squareDistance { systemParticles.squareDistancePair(positionParticle, realIndex) };
            energy += ljPairEnergy(squareDistance, particleType, systemParticles.getParticleTypeI(realIndex));
        }
    }
    return energy;
}

