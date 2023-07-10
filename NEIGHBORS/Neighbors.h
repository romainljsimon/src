/*
 * MonteCarlo.h
 *
 *  Created on: 27 oct. 2022
 *      Author: Romain Simon
 */

#ifndef NEIGHBORS_H_
#define NEIGHBORS_H_

#include <utility>
#include <iterator>
#include <fstream>
#include <vector>
#include "../INPUT/Parameter.h"
#include "util.h"
#include "../POTENTIALS/PairPotentials.h"
#include "../POTENTIALS/BondPotentials.h"
#include "../PARTICLES/Particles.h"


class Neighbors
{

private:
	std::vector<std::vector<int>> m_neighborList {};                // Neighbor list.
    const double m_rSkin {};
    const double m_squareRSkin {};                        			// Skin radius squared.
    const std::string m_neighMethod {};                   			// Neighbor list method: "verlet" for verlet neighbor list. Any other value: no neighbor list.
    int m_updateRate {-1};
    int m_numCell {};
    double m_cellLength {};
	int m_errors { 0 };                                             // Errors of the neighbor list.
	std::vector<std::vector<double>> m_interDisplacementMatrix {};  // Inter neighbor list update displacement matrix.
	const std::string m_simulationMol {};                       			// Type of system: can be either "polymer" or "atomic".
    const int m_nParticles {};
    const std::vector<double> m_maxSquareRcArray{};
    const std::vector<double> m_threshArray{};

public:
    // Polymer Neighbors. The neighbor list does NOT include bonded neighbors.
    Neighbors (param::Parameter param, const Particles& systemParticles, const PairPotentials& systemPairPotentials,
               const BondPotentials& systemBondPotentials)

            : m_simulationMol (param.get_string("simType"))
            , m_rSkin { param.get_double( "rSkin") }
            , m_squareRSkin {pow (m_rSkin, 2 ) }
            , m_neighMethod (param.get_string( "neighMethod", "verlet"))
            , m_maxSquareRcArray (initializeMaxRc( systemPairPotentials))
            , m_threshArray (initializeThreshArray(param, m_maxSquareRcArray))
            , m_numCell { static_cast<int>(systemParticles.getLengthCube() / m_rSkin) }
            , m_cellLength { systemParticles.getLengthCube() / static_cast<double>(m_numCell)}
            , m_nParticles {systemParticles.getNParticles()}

    {

        m_interDisplacementMatrix.resize(m_nParticles, std::vector<double>(3, 0));
        m_neighborList.resize(m_nParticles);
        WOWcreateNeighborList(systemParticles, systemBondPotentials);
        /***
        for(auto const& x : m_neighborList)
        {
            for(auto const& y : x)
            {
                std::cout << y << " ";
            }
            std::cout << "\n";
        }
         ***/
    }

    // Atomic constructor
    Neighbors (param::Parameter param, const Particles& systemParticles,
               const PairPotentials& systemPairPotentials)

            : m_simulationMol (param.get_string("simType"))
            , m_rSkin { param.get_double( "rSkin") }
            , m_squareRSkin {pow (m_rSkin, 2 ) }
            , m_neighMethod (param.get_string( "neighMethod", "verlet"))
            , m_maxSquareRcArray (initializeMaxRc( systemPairPotentials))
            , m_threshArray (initializeThreshArray(param, m_maxSquareRcArray))
            , m_numCell { static_cast<int>(systemParticles.getLengthCube() / m_rSkin) }
            , m_cellLength { systemParticles.getLengthCube() / static_cast<double>(m_numCell)}
            , m_nParticles {systemParticles.getNParticles()}

    {

        m_interDisplacementMatrix.resize(m_nParticles, std::vector<double>(3, 0));
        m_neighborList.resize(m_nParticles);
        createNeighborList(systemParticles);


    }


    static std::vector<double> initializeMaxRc(const PairPotentials& systemPairPotentials)
    {
        int nParticleTypes { systemPairPotentials.getParticleTypes() };
        std::vector<double> maxSquareRcArray(nParticleTypes, 0) ;

        for (int i=1; i<=nParticleTypes; i++)
        {
            for (int j=i; j<=nParticleTypes; j++)
            {
                double squareRcIJ { systemPairPotentials.getSquareRcIJ(i, j)};
                if (squareRcIJ > maxSquareRcArray[i-1])
                {
                    maxSquareRcArray[i-1] = squareRcIJ;
                }
                if (squareRcIJ > maxSquareRcArray[j-1])
                {
                    maxSquareRcArray[j-1] = squareRcIJ;
                }

            }
        }
        return maxSquareRcArray;
    };

    static std::vector<double> initializeThreshArray(param::Parameter param, const std::vector<double>& maxSquareRcArray)
    {
        double rSkin {param.get_double("rSkin")};
        std::vector<double> threshArray ;

        for (auto const &elt: maxSquareRcArray)
        {
            double maxRc { pow(elt, 1./2.) };
            const double rootThresh { (rSkin - maxRc)  / 2.};
            threshArray.push_back(rootThresh * rootThresh);
        }
        return threshArray;
    };

	void checkInterDisplacement(const Particles& systemParticles);


    [[nodiscard]] int cellTest(int indexCell) const;

    void createNeighborList(const Particles& systemParticles);
    void WOWcreateNeighborList(const Particles& systemParticles, const BondPotentials& systemBondPotentials);

    [[nodiscard]] std::vector<std::vector<std::vector<std::vector<int>>>>
    createCellList(const Particles &systemParticles) const;

    void createCellNeighbors(const Particles& systemParticles,
                                        const std::vector<std::vector<int>>& oldNeighborList,
                                        const std::vector<std::vector<std::vector<std::vector<int>>>>& cellList,
                                        int xCell, int yCell, int zCell, const bool& checkNeigh);

    void createSamePairCellNeighbor(const Particles& systemParticles,
                                    const std::vector<std::vector<int>>& oldNeighborList,
                                    const std::vector<int>& cellParticles, const bool& checkNeigh);

    void createDiffPairCellNeighbor(const Particles& systemParticles,
                                    const std::vector<std::vector<int>>& oldNeighborList,
                                    const std::vector<int>& cellParticles,
                                    const std::vector<int>& testNeighbors, const bool& checkNeigh);

    void checkNeighborError(const Particles& systemParticles, const std::vector<std::vector<int>> &oldNeighborList,
                            const double &squareDistance, const int &indexI, const int &indexJ);



    void updateIJNeighbor(const Particles& systemParticles, const std::vector<std::vector<int>>& oldNeighborList,
                          const std::vector<double>& positionParticle, const int& indexI, const int& indexJ,
                          const bool& checkNeigh);


    [[nodiscard]] const std::vector<int>& getNeighborIList(int particleIndex) const;

    void updateInterDisplacement(const int &indexTranslation, const std::vector<double> &vectorTranslation);

    void sortNeighborList();

    void createNeighborList(const Particles &systemParticles, const BondPotentials &systemBondPotentials);

    void createCellNeighbors(const Particles &systemParticles, const BondPotentials &systemBondPotentials,
                             const std::vector<std::vector<int>> &oldNeighborList,
                             const std::vector<std::vector<std::vector<std::vector<int>>>> &cellList, int xCell,
                             int yCell,
                             int zCell, const bool &checkNeigh);

    void createSamePairCellNeighbor(const Particles &systemParticles, const BondPotentials &systemBondPotentials,
                                    const std::vector<std::vector<int>> &oldNeighborList,
                                    const std::vector<int> &cellParticles, const bool &checkNeigh);

    void createDiffPairCellNeighbor(const Particles &systemParticles, const BondPotentials &systemBondPotentials,
                                    const std::vector<std::vector<int>> &oldNeighborList,
                                    const std::vector<int> &cellParticles, const std::vector<int> &testNeighbors,
                                    const bool &checkNeigh);

    void checkInterDisplacement(const Particles &systemParticles, const BondPotentials &systemBondPotentials);

    [[nodiscard]] int getUpdateRate() const;

    [[nodiscard]] int getErrors() const;
};


#endif /* NEIGHBORS_H_ */
