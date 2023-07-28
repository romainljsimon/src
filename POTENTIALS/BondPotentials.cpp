

#include <vector>
#include <cmath>
#include "BondPotentials.h"

int BondPotentials::getIndexIJ(const int& i, const int& j) const
{
    int indexI{ i };
    int indexJ{ j };

    if (i > j)
    {
        indexI = j;
        indexJ = i;
    }

    const int indexIJ {3 * ( 2 * indexJ +  indexI + 2 * m_particleTypes * ( indexI - 1) -  indexI * indexI  - 2)};
    return indexIJ;
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
    auto it {m_bondPotentials.begin() + indexIJ};
    const double& squareR0IJ {*it};
    ++it;
    double energy {0.};



    if (squareDistance >= squareR0IJ)
    {
        return std::numeric_limits<double>::infinity();
    }
    else
    {
        const double& feneKI {*it};
        ++it;

        // std::cout << feneKI << "  " << squareR0IJ << "\n";
        energy += -0.5 * feneKI * squareR0IJ * std::log(1. - squareDistance / squareR0IJ);
    }

    double rcSquareIJ {*it};
    ++it;

    if (squareDistance < rcSquareIJ)
    {

        const double& epsilonIJ {*it};
        ++it;
        const double& squareSigmaIJ {*it};
        ++it;
        const double& shiftIJ {*it};
        ++it;
        const double rapSquare { squareSigmaIJ / squareDistance };
        const double rapSix {rapSquare * rapSquare * rapSquare};
        //std::cout << epsilonIJ << "  " << squareSigmaIJ << "  " << shiftIJ << "\n";
        energy += 4. * epsilonIJ * rapSix * ( rapSix - 1.) + 4. * epsilonIJ * shiftIJ;
    }
    return energy;
}
