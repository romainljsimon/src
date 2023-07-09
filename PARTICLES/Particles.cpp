

#include <vector>
#include "Particles.h"


int Particles::getNParticles() const
{
    return m_nParticles;
}

double Particles::getLengthCube() const
{
    return m_lengthCube;
}

double Particles::getHalfLengthCube() const
{
    return m_halfLengthCube;
}

std::vector<double> Particles::getPositionI(const int& i) const
{
    return m_positionArray[i];
}
int Particles::getParticleTypeI(const int& i) const
{
    return m_particleTypeArray[i];
}

int Particles::getMoleculeTypeI(const int& i) const
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
    int typeJ {getParticleTypeI(i)};
    m_particleTypeArray[j] = m_particleTypeArray[i];
    m_particleTypeArray[i] = typeJ;

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
    return squareDistancePair(m_positionArray[indexI], m_positionArray[indexJ]);
}

double Particles::squareDistancePair(const std::vector<double>& positionI, const int& indexJ) const
{
    return squareDistancePair(positionI, m_positionArray[indexJ]);
}

double Particles::squareDistancePair(const int& indexI, const std::vector<double>& positionJ) const
{
    return squareDistancePair(m_positionArray[indexI], positionJ);
}