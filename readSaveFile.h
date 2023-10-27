/*
 * readSaveFile.h
 *
 *  Created on: 7 oct. 2022
 *      Author: Romain Simon
 */

#ifndef READSAVEFILE_H_
#define READSAVEFILE_H_


#include <vector>
#include <fstream>

void saveInXYZ(const std::vector<std::vector<double>>& positionArray, const std::vector<double>& radiusArray,
               const std::vector<int>& moleculeType, const double& lengthCube,  const std::string& path);
void saveDoubleTXT(const double& number, const std::string& path);

template <typename InputIt>
void saveDisplacement(InputIt dispItBegin1, const int& lenVec, const int& nDims, const std::string& path)
{
    /*
     * This function saves in a .xyz file the radius and the position of each particle.
      */
    std::ofstream fOut(path);
    std::string space {" "};
    for(auto it=dispItBegin1; it<dispItBegin1+lenVec; it+=nDims)
    {
        for (int i=0; i<nDims; i++)
        {
            fOut << *(it+i);
            if (i<(nDims-1))
            {
                fOut << space;
            }
        }
        fOut << "\n";
    }

    fOut.close();
}

std::vector<std::vector<int>> readBondsTXT(const std::string& path);
void saveDoubleIntTXT(const double& number1, const int& number2, const std::string& path);


#endif /* READSAVEFILE_H_ */
