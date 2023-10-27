//
// Created by rsimon on 25/10/23.
//

#ifndef BASEMOL_H
#define BASEMOL_H

#include "Molecules.h"
#include "../util.h"

class BaseMol
{
private:
    // const Molecules m_systemMolecules;
    std::vector<double> m_molBase;
    int m_lenBase;
    using BaseIterator = std::vector<double>::const_iterator;

public:
    // Monte Carlo constructor
    explicit BaseMol (const Molecules& systemMolecules)
    : m_molBase(createMolBase(systemMolecules))
    , m_lenBase(9 * systemMolecules.getNMolecules())
    {
    }

    static std::vector<double> createMolBase(const Molecules& systemMolecules)
    {
        // get unwrapped positions for molecules
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
            double normVectorY {getNormVector(vectorY.begin(), vectorY.end())};
            divVector(vectorY.begin(), 3, normVectorY);
            std::vector<double> vectorZ { crossProduct(vector2.begin(), vector1.begin(), nDims)};
            double normVectorZ {getNormVector(vectorZ.begin(), vectorZ.end())};
            divVector(vectorZ.begin(), 3, normVectorZ);
            std::vector<double> vectorX { crossProduct(vectorY.begin(), vectorZ.begin(), nDims)};

            molBase.insert(molBase.end(), vectorX.begin(), vectorX.end());
            molBase.insert(molBase.end(), vectorY.begin(), vectorY.end());
            molBase.insert(molBase.end(), vectorZ.begin(), vectorZ.end());
        }
    return molBase;
    }

    template<typename InputIt1, typename  InputIt2>
    std::vector<double> createRotationMatrix(InputIt1 molBaseBegin1, InputIt2 molBaseBegin2)
    {
        std::vector<double> rot;
        rot.reserve(9);
        for (auto it1=molBaseBegin1; it1 < molBaseBegin1+9; it1+=3)
        {
            for (auto it2=molBaseBegin2; it2 < molBaseBegin2+9; it2+=3)
            {
                double elt = (*it1) * (*it2) + (*(it1+1)) * (*(it2+1)) + (*(it1+2)) * (*(it2+2));
                rot.push_back(elt);
            }
        }
        return rot;
    }

    template<typename InputIt1, typename  InputIt2>
    std::vector<double> createAxisAngle(InputIt1 molBaseBegin1, InputIt2 molBaseBegin2)
    {
        // find axis direction and normalize vector
        std::vector<double> axisAngle(3);
        std::vector<double> rotMatrix {createRotationMatrix(molBaseBegin1, molBaseBegin2)};
        axisAngle[0] = rotMatrix[7] - rotMatrix[5];
        axisAngle[1] = rotMatrix[2] - rotMatrix[6];
        axisAngle[2] = rotMatrix[3] - rotMatrix[1];

        double normAxisAngle {getNormVector(axisAngle.begin(), axisAngle.end())};

        if (normAxisAngle>0)
        {
            divVector(axisAngle.begin(), 3, normAxisAngle);
            // find angle
            std::vector<double> kn {0, -axisAngle[2], axisAngle[1], axisAngle[2], 0, -axisAngle[0], -axisAngle[1], axisAngle[0], 0};
            std::vector<double> matMul {matrixMultiplication(kn.begin(), rotMatrix.begin(), 3, 3, 3)};
            double cosAngle {traceMatrix(rotMatrix.begin(), 3) / 2. - 1./2.};
            double sinAngle { -traceMatrix(matMul.begin(), 3) / 2.};
            double sign {1};
            if (sinAngle < 0)
            {
                //std::cout << sinAngle*sinAngle + cosAngle*cosAngle << "\n";
                sign = -1;
            }
            double angle {acos(cosAngle) * sign};
            // compute axisAngle
            multVector(axisAngle.begin(), 3, angle);
        }

        return axisAngle;
    }

    template<typename InputIt>
    std::vector<double> createAngularVector(InputIt TotalMolBaseItT)
    {
        //auto molBaseTItBegin {MolBaseT.getBegin()}
        //double lenBaseMol{9 * mol};
        std::vector<double> angularVectorDt;
        angularVectorDt.reserve(m_lenBase/3);
        for (int i=0; i < m_lenBase; i+=9)
        {
            std::vector<double> axisAngle {createAxisAngle(m_molBase.begin() + i, TotalMolBaseItT + i)};
            angularVectorDt.insert(angularVectorDt.end(), axisAngle.begin(), axisAngle.end());
        }
        return angularVectorDt;
    }

    BaseIterator BaseMolBeginIt()
    {
        return m_molBase.begin();
    }

};


#endif //BASEMOL_H
