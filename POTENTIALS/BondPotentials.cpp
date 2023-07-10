

#include <vector>
#include <cmath>
#include "BondPotentials.h"

const std::vector<int>& BondPotentials::getBondsI(const int& i) const
{
    return m_bondsArray[i];
}

int BondPotentials::getIndexIJ(const int& i, const int& j) const
{
    int indexI{ i };
    int indexJ{ j };

    if (i > j)
    {
        indexI = j;
        indexJ = i;
    }

    const int indexIJ {indexJ - indexI + m_particleTypes * (indexI - 1) - ((indexI - 2) * (indexI-1)) / 2};
    return indexIJ * 6;
}

/***
std::vector<double> BondPotentials::getPotentialsIJ(const int& i, const int& j) const
{
    int indexIJ {getIndexIJ(i, j)};
    return m_bondPotentials[indexIJ];
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
double BondPotentials::feneBondEnergyIJ(const double& squareDistance, const int& particleTypeI,
                                        const int& particleTypeJ) const
{
    const int indexIJ {getIndexIJ(particleTypeI, particleTypeJ)};
    const double squareR0IJ {m_bondPotentials[indexIJ + 1]};
    double energy {0.};



    if (squareDistance >= squareR0IJ)
    {
        return std::numeric_limits<double>::infinity();
    }
    else
    {

        const double feneKI {m_bondPotentials[indexIJ]};
        // std::cout << feneKI << "  " << squareR0IJ << "\n";
        energy += -0.5 * feneKI * squareR0IJ * std::log(1. - squareDistance / squareR0IJ);
    }

    double rcSquareIJ {m_bondPotentials[indexIJ + 4]};

    if (squareDistance < rcSquareIJ)
    {

        const double epsilonIJ {m_bondPotentials[indexIJ + 2]};
        const double squareSigmaIJ {m_bondPotentials[indexIJ + 3]};
        const double shiftIJ {m_bondPotentials[indexIJ + 5]};
        const double rapSquare { squareSigmaIJ / squareDistance };
        const double rapSix {rapSquare * rapSquare * rapSquare};
        //std::cout << epsilonIJ << "  " << squareSigmaIJ << "  " << shiftIJ << "\n";
        energy += 4. * epsilonIJ * rapSix * ( rapSix - 1.) + 4. * epsilonIJ * shiftIJ;
    }
    return energy;
}

double BondPotentials::feneBondEnergyI(const int& indexParticle, const std::vector<double>& positionParticle,
                                       const Particles& systemParticles, const int& indexSkip) const
{
    double energy { 0. };

    std::vector<int> bondsI { m_bondsArray[indexParticle] };
    const int bondsISize { static_cast<int>(bondsI.size()) };


    for (int i = 0; i < bondsISize; i++)
    {
        const int realIndex = bondsI[i];

        if (realIndex != indexSkip)
        {
            //if ((realIndex == indexParticle) || (realIndex == -1))
            //{
            //    continue;
            //}
            const double squareDistance {systemParticles.squareDistancePair(positionParticle, realIndex)};
            energy += feneBondEnergyIJ(squareDistance, systemParticles.getParticleTypeI(indexParticle),
                                       systemParticles.getParticleTypeI(realIndex));
        }
    }
    return energy;
}

