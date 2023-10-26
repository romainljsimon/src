//
// Created by rsimon on 25/10/23.
//

#ifndef ROTATIONMOL_H
#define ROTATIONMOL_H

#include "../MOLECULES/Molecules.h"
#include "../util.h"

class RotationMol
{
private:
    const Molecules m_systemMolecules;
    const std::vector<double> m_molBase;

public:
    // Monte Carlo constructor
    explicit RotationMol (const Molecules& systemMolecules)
    : m_systemMolecules(systemMolecules)
    , m_molBase(createMolBase(systemMolecules))
    {
    }

    static std::vector<double> createMolBase(const Molecules& systemMolecules)
    {
        // get unwrapped positions for molecules
        //std::vector<double> positionCMMolecule {vectorfromCM(positionArray, lengthChain)};
        const int nParticles {systemMolecules.getNParticles()};
        const int nDims {systemMolecules.getNDims()};
        const double lengthCube {systemMolecules.getLengthCube()};
        const double lengthChain {3.};
        std::vector<double> molBase;
        molBase.reserve(nDims * nParticles);
        for (int indexMolecule=0; indexMolecule<nParticles; indexMolecule+=3)
        {
            // create molecular base
            std::vector<int> orderVector{systemMolecules.getOrderVector(indexMolecule)};
            std::vector<double> vectorStart{systemMolecules.getPositionI(indexMolecule + orderVector[0])};
            std::vector<double> vectorMiddle{systemMolecules.getPositionI(indexMolecule + orderVector[1])};
            std::vector<double> vectorEnd{systemMolecules.getPositionI(indexMolecule + orderVector[2])};
            std::vector<double> vector1{vectorDiff(vectorMiddle, vectorStart, nDims)};
            std::vector<double> vector2{vectorDiff(vectorEnd, vectorStart, nDims)};

            // Vector 1 is the vector linking Monomer 1 to Monomer 2
            periodicVector(vector1.begin(), nDims, lengthCube);

            // Vector 2 is the vector linking Monomer 1 to Monomer 3
            periodicVector(vector2.begin(), nDims, lengthCube);

            // VectorCM is the center of mass vector vectorCM = (vector1 + vector2)/3
            std::vector<double> vectorCM {vector1};
            sumVector(vectorCM.begin(), vector2.begin(), nDims);
            divVector(vectorCM.begin(), nDims, lengthChain);
            std::vector<double> vectorY { vectorDiff(vector1, vectorCM, nDims)};
            std::vector<double> vectorZ { crossProduct(vector2.begin(), vector1.begin(), nDims)};
            std::vector<double> vectorX { crossProduct(vectorY.begin(), vectorZ.begin(), nDims)};
            //NORMALIZE !!!!!!
            molBase.insert(molBase.end(), vectorX.begin(), vectorX.end());
            molBase.insert(molBase.end(), vectorY.begin(), vectorY.end());
            molBase.insert(molBase.end(), vectorZ.begin(), vectorZ.end());
        }
    return molBase;
    }
};


#endif //ROTATIONMOL_H
