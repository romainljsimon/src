//
// Created by rsimon on 25/10/23.
//

#include "RotationMol.h"


std::vector<double> createRotationMatrix(const std::vector<double>& molBase1, const std::vector<double>& molBase2)
{
    std::vector<double> rot;
    rot.reserve(9);
    for (auto it1=molBase1.begin(); it1 != molBase1.end(); it1+=3)
    {
        for (auto it2=molBase2.begin(); it2 != molBase2.end(); it2+=3)
        {
            double elt = (*it1) * (*it2) + (*(it1+1)) * (*(it2+1)) + (*(it1+2)) * (*(it2+2));
            rot.push_back(elt);
        }
    }
    return rot;
}

std::vector<double> createAxisAngle(const std::vector<double>& molBase1, const std::vector<double>& molBase2)
{
    std::vector<double> axisAngle(3);
    std::vector<double> rotMatrix {createRotationMatrix(molBase1, molBase2)};
    axisAngle[0] = rotMatrix[7] - rotMatrix[5];
    axisAngle[1] = rotMatrix[2] - rotMatrix[6];
    axisAngle[2] = rotMatrix[3] - rotMatrix[1];
    std::vector<double> kn {0, -axisAngle[2], axisAngle[1], axisAngle[2], 0, axisAngle[0], -axisAngle[1], axisAngle[0], 0};
    
    double cosAngle {traceMatrix(rotMatrix.begin(), 3) / 2. - 1./2.};

    return axisAngle;
}