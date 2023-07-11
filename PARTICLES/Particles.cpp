

#include <vector>
#include "Particles.h"


const int& Particles::getNParticles() const
{
    return m_nParticles;
}

const double& Particles::getLengthCube() const
{
    return m_lengthCube;
}

const double& Particles::getHalfLengthCube() const
{
    return m_halfLengthCube;
}

const std::vector<double>& Particles::getPositionI(const int& i) const
{
    return m_positionArray[i];
}
const int& Particles::getParticleTypeI(const int& i) const
{
    return m_particleTypeArray[i];
}

const int& Particles::getMoleculeTypeI(const int& i) const
{
    return m_moleculeTypeArray[i];
}

void Particles::updatePositionI(const int& i, const std::vector<double>& newPosition)
{
    m_positionArray[i] = newPosition;
}

void Particles::swapParticleTypesIJ(const int& i, const int& j, const int& typeI, const int& typeJ)
{
    m_particleTypeArray[i] = typeJ;
    m_particleTypeArray[j] = typeI;
}

void Particles::swapParticleTypesIJ(const int& i, const int& j)
{
    //std::cout << i << "  " << j << "\n";
    //std::cout << m_particleTypeArray[i] << "  " << m_particleTypeArray[j] << "\n";

    const int typeJ {getParticleTypeI(j)};
    //std::cout << m_particleTypeArray[i] << m_particleTypeArray[j] << "\n";
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
double Particles::squareDistancePair(const std::vector<double>& positionI, const std::vector<double>& positionJ) const
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

double Particles::squareDistancePair(const int& indexI, const int& indexJ) const
{
    const std::vector<double>& positionI {m_positionArray[indexI]};
    const std::vector<double>& positionJ {m_positionArray[indexJ]};
    return squareDistancePair(positionI, positionJ);
}

double Particles::squareDistancePair(const std::vector<double>& positionI, const int& indexJ) const
{
    const std::vector<double>& positionJ {m_positionArray[indexJ]};
    return squareDistancePair(positionI, positionJ);
}

double Particles::squareDistancePair(const int& indexI, const std::vector<double>& positionJ) const
{
    const std::vector<double>& positionI {m_positionArray[indexI]};
    return squareDistancePair(positionI, positionJ);
}


std::vector<double> Particles::periodicBC(std::vector<double> positionParticle)
/*
 *This function is an implementation of the periodic Boundary conditions.
 *If a particle gets out of the simulation box from one of the sides, it gets back in the box from the opposite side.
 */
{
    m_newFlags = {0, 0, 0};

    for (int i = 0; i < 3; i++)
    {
        double& positionI {positionParticle[i]};
        if (positionI < 0)
        {
            positionI += m_lengthCube;
            m_newFlags[i] -= 1;
        }
        else if (positionI > m_lengthCube)
            positionI -= m_lengthCube;
            m_newFlags[i] += 1;
    }
    return positionParticle;
}

void Particles::updateFlags(const int& indexParticle)
{
    for (int i = 0; i < 3; i++)
    {
        m_flagsArray[indexParticle + i] += m_newFlags[i];
    }
}


void Particles::saveInXYZ(const std::string& path) const
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
    fOut << " Properties=type:I:1:radius:R:1:pos:R:3";
    fOut << "\n";
    int it {0};
    for(auto const& x : m_positionArray)
    {

        int j { 0 };
        for(auto const& i : x)
        {

            if (j == 2)
            {
                fOut << i;
                fOut << " ";
                fOut << m_flagsArray[3*it];
                fOut << " ";
                fOut << m_flagsArray[3*it+1];
                fOut << " ";
                fOut << m_flagsArray[3*it+2];
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
