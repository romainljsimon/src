

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
    const int indexIJ {2 * ( 2 * indexJ +  indexI + 2 * m_nParticleTypes * ( indexI - 1) -  indexI * indexI  - 2)};
    return indexIJ;
}

int PairPotentials::getParticleTypes() const
{
    return m_nParticleTypes;
}

double PairPotentials::getSquareRcIJ(const int& i, const int& j) const
{
    const int rcIndex {getIndexIJ(i, j)};
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
    auto it {m_pairPotentials.begin() + indexIJ};
    const double& rcSquareIJ { *it};
    ++it;

    if (squareDistance > rcSquareIJ)
    {

        return 0.;
    }

    else
    {

        const double& fourEpsilonIJ {*it};
        ++it;
        const double& squareSigmaIJ {*it};
        ++it;
        const double& shiftIJ { *it};
        ++it;
        const double rapSquare { squareSigmaIJ / squareDistance };
        const double rapSix { rapSquare * rapSquare * rapSquare};

        return fourEpsilonIJ * rapSix * ( rapSix - 1.) + fourEpsilonIJ * shiftIJ;
    }
}

