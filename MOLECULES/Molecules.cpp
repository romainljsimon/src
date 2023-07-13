

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

const std::vector<double>& Molecules::getPositionI(const int& i) const
{
    return m_positionArray[i];
}
const int& Molecules::getParticleTypeI(const int& i) const
{
    return m_particleTypeArray[i];
}

const int& Molecules::getMoleculeTypeI(const int& i) const
{
    return m_moleculeTypeArray[i];
}

void Molecules::updatePositionI(const int& i, const std::vector<double>& newPosition)
{
    m_positionArray[i] = newPosition;
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
double Molecules::squareDistancePair(const std::vector<double>& positionI, const std::vector<double>& positionJ) const
{
    double squareDistance { 0. };
    for (int i = 0; i < 3; i++)
    {
        double diff { positionI[i] - positionJ[i] };

        if (diff > m_halfLengthCube)
        {
            diff -= m_lengthCube;
        }
        else if (diff < - m_halfLengthCube)
        {
            diff += m_lengthCube;
        }

        squareDistance += diff * diff;
    }

    return squareDistance;
}

double Molecules::squareDistancePair(const int& indexI, const int& indexJ) const
{
    const std::vector<double>& positionI {m_positionArray[indexI]};
    const std::vector<double>& positionJ {m_positionArray[indexJ]};
    return squareDistancePair(positionI, positionJ);
}

double Molecules::squareDistancePair(const std::vector<double>& positionI, const int& indexJ) const
{
    const std::vector<double>& positionJ {m_positionArray[indexJ]};
    return squareDistancePair(positionI, positionJ);
}

double Molecules::squareDistancePair(const int& indexI, const std::vector<double>& positionJ) const
{
    const std::vector<double>& positionI {m_positionArray[indexI]};
    return squareDistancePair(positionI, positionJ);
}


std::vector<double> Molecules::periodicBC(std::vector<double> positionParticle)
/*
 *This function is an implementation of the periodic Boundary conditions.
 *If a particle gets out of the simulation box from one of the sides, it gets back in the box from the opposite side.
 */
{
    for (int i = 0; i < 3; i++)
    {
        //std::cout <<  m_newFlags[i] << "\n";
        const double& positionI { positionParticle[i] };

        if (positionI < 0)
        {

            positionParticle[i] += m_lengthCube;
            m_newFlags.push_back(-1);
        }
        else if (positionI > m_lengthCube)
        {
            positionParticle[i] -= m_lengthCube;
            m_newFlags.push_back(1);
        }
        else
        {
            m_newFlags.push_back(0);
        }
    }

    return positionParticle;
}

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

double Molecules::energyPairParticle(const int& indexParticle, const std::vector<double>& positionParticle,
                                     const std::vector<int>& neighborIList, const int& indexSkip) const
{
    double energy { 0. };
    const int& particleType {m_particleTypeArray[indexParticle]};


    for (int const &indexJ: neighborIList)
    {
        if (indexJ != indexSkip)
        {
            //if (realIndex == indexParticle)
            //{
            //    continue;
            //}
            const int& typeJ {m_particleTypeArray[indexJ]};
            const double squareDistance { squareDistancePair(positionParticle, indexJ) };
            energy += m_systemPairPotentials.ljPairEnergy(squareDistance, particleType, typeJ);

        }
    }
    return energy;
}

double Molecules::energyPairParticle(const int& indexParticle, const std::vector<double>& positionParticle,
                                     const std::vector<int>& neighborIList) const
{
    return energyPairParticle(indexParticle, positionParticle, neighborIList, -1);
}

double Molecules::energyPairParticle(const int& indexParticle, const std::vector<int>& neighborIList,
                                     const int& indexSkip) const
{
    const std::vector<double>& positionParticle {m_positionArray[indexParticle]};
    return energyPairParticle(indexParticle, positionParticle, neighborIList, indexSkip);
}

double Molecules::energyPairParticle(const int& indexParticle, const std::vector<int>& neighborIList) const
{
    return energyPairParticle(indexParticle, neighborIList, -1);
}

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
            const double squareDistance { squareDistancePair(positionParticle, indexJ) };
            energy += m_systemPairPotentials.ljPairEnergy(squareDistance, particleType, typeJ);

        }
    }
    return energy;
}


double Molecules::energyPairParticleExtraMolecule(const int& indexParticle, const std::vector<int>& neighborIList,
                                     const int& typeMoleculeI) const
{
    const std::vector<double>& positionParticle {m_positionArray[indexParticle]};
    return energyPairParticleExtraMolecule(indexParticle, positionParticle, neighborIList, typeMoleculeI);
}


double Molecules::feneBondEnergyI(const int& indexParticle, const std::vector<double>& positionParticle,
                                  const int& indexSkip) const
{
    double energy { 0. };

    const std::vector<int>& bondsI { m_bondsArray[indexParticle] };
    const int& particleTypeI { m_particleTypeArray[indexParticle] };
    for (int const &indexJ: bondsI)
    {

        if (indexJ != indexSkip)
        {
            const double squareDistance { squareDistancePair(positionParticle, indexJ)};
            const int& particleTypeJ { m_particleTypeArray[indexJ] };
            energy += m_systemBondPotentials.feneBondEnergyIJ(squareDistance, particleTypeI, particleTypeJ);
        }
    }
    return energy;
}

double Molecules::feneBondEnergyI(const int& indexParticle, const int& indexSkip) const
{
    const std::vector<double>& positionParticle {m_positionArray[indexParticle]};
    return feneBondEnergyI(indexParticle, positionParticle, indexSkip);
}

double Molecules::feneBondEnergyI(const int& indexParticle) const
{
    return feneBondEnergyI(indexParticle, -1);
}


double Molecules::energyParticleMolecule(const int& indexParticle, const std::vector<double>& positionParticle,
                                 const std::vector<int>& neighborIList, const int& indexSkip) const
{
    double energy { 0. };
    energy += energyPairParticle(indexParticle, positionParticle, neighborIList, indexSkip);
    energy += feneBondEnergyI(indexParticle, positionParticle, indexSkip);
    return energy;
}

double Molecules::energyParticleMolecule(const int& indexParticle, const std::vector<double>& positionParticle,
                                         const std::vector<int>& neighborIList) const
{
    return energyParticleMolecule(indexParticle, positionParticle, neighborIList, -1);
}

double Molecules::energyParticleMolecule(const int& indexParticle, const std::vector<int>& neighborIList,
                                         const int& indexSkip) const
{

    const std::vector<double>& positionParticle {m_positionArray[indexParticle]};

    return energyParticleMolecule(indexParticle, positionParticle, neighborIList, indexSkip);
}


double Molecules::energyParticleMolecule(const int& indexParticle, const std::vector<int>& neighborIList) const
{
    return energyParticleMolecule(indexParticle, neighborIList, -1);
}


double Molecules::energySystemMolecule(const Neighbors& systemNeighbors) const
{
    double energy { 0. };

    for (int indexParticle = 0; indexParticle < m_nParticles; indexParticle++) //Outer loop for rows
    {
        const std::vector<int> &neighborIList{systemNeighbors.m_neighborList[indexParticle]};
        energy += energyParticleMolecule(indexParticle, neighborIList) / 2.;
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
    fOut << m_nParticles;
    fOut << "\n";
    std::string lengthStr = std::to_string(m_lengthCube);
    fOut << "Lattice=";
    fOut << '"';
    fOut << lengthStr;
    fOut << " 0.0 0.0 0.0 ";
    fOut << lengthStr;
    fOut << " 0.0 0.0 0.0 ";
    fOut << lengthStr;
    fOut << '"';
    fOut << " Properties=molecule_type:S:1:type:I:1:pos:R:3:";
    fOut << "\n";
    int it {0};
    for(auto const& x : m_positionArray)
    {

        int j { 0 };
        for(auto const& i : x)
        {

            if (j == 2)
            {
                const int flagIt {3*it};
                fOut << i;
                fOut << " ";
                fOut << m_flagsArray[flagIt];
                fOut << " ";
                fOut << m_flagsArray[flagIt+1];
                fOut << " ";
                fOut << m_flagsArray[flagIt+2];
                fOut << "\n";
            }
            else if (j==0)
            {
                fOut << m_moleculeTypeArray[it];
                fOut << " ";
                fOut << m_particleTypeArray[it];
                fOut << " ";
                fOut << i;
                fOut << " ";
            }
            else
            {
                fOut << i;
                fOut << " ";
            }

            j += 1;

        }
        it += 1;
    }
    fOut.close();
}

