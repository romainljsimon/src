

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
    const double& rcSquareIJ { m_pairPotentials[indexIJ + 2] };

    if (squareDistance > rcSquareIJ)
    {

        return 0.;
    }

    else
    {

        const double& epsilonIJ { m_pairPotentials[indexIJ + 0]};

        const double& squareSigmaIJ { m_pairPotentials[indexIJ + 1]};
        const double& shiftIJ { m_pairPotentials[indexIJ + 3] };
        const double rapSquare { squareSigmaIJ / squareDistance };
        const double rapSix { rapSquare * rapSquare * rapSquare};

        return 4. * epsilonIJ * rapSix * ( rapSix - 1.) + 4. * epsilonIJ * shiftIJ;
    }
}

