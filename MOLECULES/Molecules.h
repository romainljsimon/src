#ifndef MOLECULES_H_
#define MOLECULES_H_

#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <vector>
#include "../INPUT/Parameter.h"
#include "../POTENTIALS/PairPotentials.h"
#include "../POTENTIALS/BondPotentials.h"

class Neighbors;

class Molecules
{
    friend class Neighbors;

public:
    // POTENTIALS constructor

    const PairPotentials m_systemPairPotentials {};
    const BondPotentials m_systemBondPotentials {};
    const std::vector<std::vector<int>> m_bondsArray {};
    const int m_nParticles {};
    const double m_lengthCube {};
    const double m_halfLengthCube {};
    std::vector<int> m_newFlags;
    std::vector<int> m_flagsArray {};
    std::vector<std::vector<double>> m_positionArray {};
    std::vector<int> m_particleTypeArray {};
    std::vector<int> m_moleculeTypeArray {};

    Molecules (param::Parameter param, PairPotentials  systemPairPotentials,
               BondPotentials  systemBondPotentials, const std::string& path)

        : m_nParticles ( initializeNumber(path))
        , m_systemPairPotentials(std::move(systemPairPotentials))
        , m_systemBondPotentials(std::move(systemBondPotentials))
        , m_lengthCube {pow (static_cast<double>(m_nParticles) / param.get_double("density"), 1. / 3.) }
        , m_halfLengthCube (0.5 * m_lengthCube)
        , m_bondsArray(initializeBondsArray())
    {
        m_flagsArray.resize(3 * m_nParticles, 0);
        initializeParticles(path);
    }

    static std::vector<std::vector<int>> initializeBondsArray()
    {
        std::ifstream infile("./bonds.txt");

        if (!infile.is_open())
            std::cout << "Error opening file";


        int nParticles{};
        infile >> nParticles;

        int nBonds{};
        infile >> nBonds;

        std::vector<std::vector<int>> bondsArray(nParticles);

        for (int r = 0; r < nBonds; r++) //Outer loop for rows
        {

            int indexI {};
            infile >> indexI;

            int indexJ {};
            infile >> indexJ;

            bondsArray[indexI].push_back(indexJ);
            bondsArray[indexJ].push_back(indexI);

        }
        infile.close();
        /***

        std::cout << bondsArray.size() << "\n";
         ***/
        return bondsArray;
    }

    static int initializeNumber(const std::string& path)
    {
        std::ifstream infile(path);

        if (!infile.is_open())
            std::cout << "Error opening file" ;

        int nParticles {};
        infile >> nParticles;
        infile.close();
        return nParticles;
    }

    void initializeParticles(const std::string& path)
    {
        std::ifstream infile(path);

        if (!infile.is_open())
            std::cout << "Error opening file" ;

        int row{};
        infile >> row;

        int col { 5 };

        std::string line;
        getline(infile, line);
        getline(infile, line);

        m_moleculeTypeArray.resize(row);
        m_particleTypeArray.resize(row);
        m_positionArray.resize(row, std::vector<double>(3));
        //std::vector<std::vector <double>> positionArray(row, std::vector<double>(3));
        //std::vector<double> typeArray (row);
        //std::vector<int> moleculeType (row , 1);

        //Defining the loop for getting INPUT from the file

        for (int r = 0; r < row; r++) //Outer loop for rows
        {

            for (int c = 0; c < col; c++) //inner loop for columns
            {

                if (c == 0)
                {
                    infile >> m_moleculeTypeArray[r];
                }
                else if (c == 1)
                {

                    infile >> m_particleTypeArray[r];
                }

                else
                {
                    infile >> m_positionArray[r][c - (col - 3)]; //Take INPUT from file and put into positionArray
                }
            }
        }
        infile.close();
    }

    void updateFlags(const int& indexParticle);

    [[nodiscard]] const int& getNParticles() const;

    [[nodiscard]]  const std::vector<double>& getPositionI(const int& i) const;

    [[nodiscard]]  const int& getParticleTypeI(const int &i) const;

    [[nodiscard]] const int& getMoleculeTypeI(const int &i) const;

    void swapParticleTypesIJ(const int &i, const int &j, const int &typeI, const int &typeJ);

    void updatePositionI(const int &i, const std::vector<double> &newPosition);

    [[nodiscard]] double squareDistancePair(const std::vector<double> &positionA,
                                            const std::vector<double> &positionB) const;

    [[nodiscard]] double squareDistancePair(const int &indexI, const int &indexJ) const;

    [[nodiscard]] double squareDistancePair(const std::vector<double> &positionI, const int &indexJ) const;

    [[nodiscard]] double squareDistancePair(const int &indexI, const std::vector<double> &positionJ) const;

    [[nodiscard]] const double& getLengthCube() const;

    [[nodiscard]] const double& getHalfLengthCube() const;

    void saveInXYZ(const std::string& path) const;

    void swapParticleTypesIJ(const int &i, const int &j);

    std::vector<double> periodicBC(std::vector<double> positionParticle);

    [[nodiscard]] double energyPairParticle(const int& indexParticle, const std::vector<double>& positionParticle,
                                            const std::vector<int>& neighborIList, const int& indexSkip) const;

    [[nodiscard]] double energyPairParticle(const int& indexParticle, const std::vector<double>& positionParticle,
                                            const std::vector<int>& neighborIList) const;

    [[nodiscard]] double energyPairParticle(const int& indexParticle, const std::vector<int>& neighborIList,
                                            const int& indexSkip) const;

    [[nodiscard]] double energyPairParticle(const int& indexParticle, const std::vector<int>& neighborIList) const;

    [[nodiscard]] double energySystemMolecule(const Neighbors &systemNeighbors) const;
    [[nodiscard]] double feneBondEnergyI(const int& indexParticle, const std::vector<double>& positionParticle,
                                         const int& indexSkip) const;

    [[nodiscard]] double feneBondEnergyI(const int &indexParticle) const;

    [[nodiscard]] double feneBondEnergyI(const int &indexParticle, const int &indexSkip) const;

    [[nodiscard]] double energyParticleMolecule(const int &indexParticle, const std::vector<int> &neighborIList) const;

    [[nodiscard]] double energyParticleMolecule(const int &indexParticle, const std::vector<double> &positionParticle,
                          const std::vector<int> &neighborIList, const int &indexSkip) const;

    [[nodiscard]] double
    energyParticleMolecule(const int &indexParticle, const std::vector<int> &neighborIList, const int &indexSkip) const;

    [[nodiscard]] double energyParticleMolecule(const int &indexParticle, const std::vector<double> &positionParticle,
                                  const std::vector<int> &neighborIList) const;

    void reinitializeFlags();

    [[nodiscard]] double energyPairParticleExtraMolecule(const int& indexParticle, const std::vector<double>& positionParticle,
                                                      const std::vector<int>& neighborIList, const int& typeMoleculeI) const;

    [[nodiscard]] double energyPairParticleExtraMolecule(const int& indexParticle, const std::vector<int>& neighborIList, const int& typeMoleculeI) const;
};

#endif /* MOLECULES_H_ */
