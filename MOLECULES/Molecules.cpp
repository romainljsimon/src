

#include <vector>
#include "Molecules.h"
#include "../NEIGHBORS/Neighbors.h"


const int& Molecules::getNParticles() const
{
    return m_nParticles;
}


const double& Molecules::getLengthCube() const
{
    return m_lengthCube;
}

const double& Molecules::getHalfLengthCube() const
{
    return m_halfLengthCube;
}

const int& Molecules::getNDims() const
{
    return m_nDims;
}


std::vector<double> Molecules::getPositionI(const int& i) const
{
    std::vector<double> positionI (getPosItBeginI(i),
                                   getPosItEndI(i));
    return positionI;
}



const int& Molecules::getParticleTypeI(const int& i) const
{
    return m_particleTypeArray[i];
}

const int& Molecules::getMoleculeTypeI(const int& i) const
{
    return m_moleculeTypeArray[i];
}


void Molecules::swapParticleTypesIJ(const int& i, const int& j, const int& typeI, const int& typeJ)
{
    m_particleTypeArray[i] = typeJ;
    m_particleTypeArray[j] = typeI;
}

void Molecules::swapParticleTypesIJ(const int& i, const int& j)
{
    //std::cout << i << "  " << j << "\n";
    //std::cout << m_particleTypeArray[i] << "  " << m_particleTypeArray[j] << "\n";

    const int typeJ {getParticleTypeI(j)};
    m_particleTypeArray[j] = m_particleTypeArray[i];
    //std::cout << m_particleTypeArray[i] << m_particleTypeArray[j] << "\n";
    m_particleTypeArray[i] = typeJ;
    //std::cout << m_particleTypeArray[i] << m_particleTypeArray[j] << "\n";
}


/*******************************************************************************
 * This function calculates the distance between two particles considering the
 * periodic boundary conditions.
 *
 * @param positionA particle A's position.
 *        particleB particle B's position.
 *
 * @return distance between two particles.
 ******************************************************************************/


/***
double Molecules::squareDistancePair(const int& indexI, const int& indexJ) const
{
    const std::vector<double>& positionI {getPositionI(indexI)};
    const std::vector<double>& positionJ { getPositionI(indexJ)};
    return squareDistancePair(positionI, positionJ);
}

double Molecules::squareDistancePair(auto positionI, const int& indexJ) const
{
    const std::vector<double> positionJ {getPositionI(indexJ)};
    return squareDistancePair(positionI, positionJ);
}

double Molecules::squareDistancePair(const int& indexI, const std::vector<double>& positionJ) const
{
    const std::vector<double>& positionI {getPositionI(indexI)};
    return squareDistancePair(positionI, positionJ);
}
***/



void Molecules::updateFlags(const int& indexParticle)
{
    int indexFlag {3 * indexParticle};
    for(auto const& flag : m_newFlags)
    {

        m_flagsArray[indexFlag] += flag;
        indexFlag++;
    }
}

void Molecules::reinitializeFlags()
{
    m_newFlags.clear();
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



double Molecules::energyPairParticleExtraMolecule(const int& indexParticle, const std::vector<double>& positionParticle,
                                     const std::vector<int>& neighborIList, const int& typeMoleculeI) const
{
    double energy { 0. };
    const int& particleType {m_particleTypeArray[indexParticle]};

    for (int const &indexJ: neighborIList)
    {
        const int& typeMoleculeJ {m_moleculeTypeArray[indexJ]};

        if (typeMoleculeJ != typeMoleculeI)
        {
            //if (realIndex == indexParticle)
            //{
            //    continue;
            //}
            const int& typeJ {m_particleTypeArray[indexJ]};
            const double squareDistance { squareDistancePair(positionParticle.begin(), m_positionArray.begin() + m_nDims * indexJ)  };
            energy += m_systemPairPotentials.ljPairEnergy(squareDistance, particleType, typeJ);

        }
    }
    return energy;
}


double Molecules::energyPairParticleExtraMolecule(const int& indexParticle, const std::vector<int>& neighborIList,
                                     const int& typeMoleculeI) const
{
    const std::vector<double>& positionParticle {getPositionI(indexParticle)};
    return energyPairParticleExtraMolecule(indexParticle, positionParticle, neighborIList, typeMoleculeI);
}


double Molecules::energySystemMolecule(const Neighbors& systemNeighbors) const
{
    double energy { 0. };

    for (int indexParticle = 0; indexParticle < m_nParticles; indexParticle++) //Outer loop for rows
    {
        const auto& neighItBegin { systemNeighbors.getNeighItBeginI(indexParticle) };
        const int& lenNeigh { systemNeighbors.getLenIndexBegin(indexParticle)};
        const double particleEnergy { energyParticleMolecule(indexParticle, neighItBegin, lenNeigh) };
        energy += particleEnergy / 2.;

     }
    return energy;
}

/*******************************************************************************
 * This function calculates the total potential energy of the system considering
 * that particles are Lennard-Jones particles.
 *
 * @param positionArray array of particle positions.
 *        typeArray array of particle diameter.
 *        squareRc square of the cut off radius.
 *        lengthCube simulation box length.
 *
 * @return System's total energy.
 ******************************************************************************/
/***
double Particles::energySystem(const Neighbors& systemNeighbors) const
{
    double energy { 0. };

    for (int indexParticle = 0; indexParticle < m_nParticles; indexParticle++) //Outer loop for rows
    {
        const std::vector<int>& neighborIList {systemNeighbors.m_neighborList[indexParticle]};
        energy += energyPairParticle(indexParticle, neighborIList) / 2.;
    }
    return energy;
}
***/

void Molecules::saveInXYZ(const std::string& path) const
{

    /*
     * This function saves in a .xyz file the radius and the position of each particle.
      */


    std::ofstream fOut(path);

    fOut << m_saveHeaderString;
    int it {0};
    std::string space {" "};
    for(auto const& x : m_positionArray)
    {

        if ((it%m_nDims)==0)
        {
            fOut << m_moleculeTypeArray[it];
            fOut << space;
            fOut << m_particleTypeArray[it];
            fOut << space;

        }

        fOut << x;
        fOut << space;

        if ((it%m_nDims)==m_nDims-1)
        {
            fOut << m_flagsArray[it];
            fOut << space;
            fOut << m_flagsArray[it+1];
            fOut << space;
            fOut << m_flagsArray[it+2];
            fOut << "\n";
        }
        it++;

    }
    fOut.close();
}
