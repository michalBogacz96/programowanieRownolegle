#ifndef SIMULATION_H

#define SIMULATION_H

#include<iostream>
#include"EnergyCalculator.h"
#include"RandomNumberGenerator.h"
#include"MonteCarlo.h"
#include <omp.h>
#include <vector>

using namespace std;

class Simulation {
        private:
                struct drand48_data* structure;

                MonteCarlo *mc;
                RandomNumberGenerator *rng;
                EnergyCalculator *energyCalculator;
                int size;
                int reducedSize;
                int dataToChange;
                int maxNeighbours;
                double maxChange;
                double *delta;
                int *rows;
                int *cols;
                double *data;
                double Etot;


                void changeData();
                void changeDataUndo();
                void generateDataChange();

                int similarNeighbours( int col, int row, int delta, double limit );
	public:
                struct drand48_data state;
                int stepNumber;
                int maxNumberOfStepsForTheSameRandomData;
                int numberOfDataToGenerate;
                std::vector<double> rowResults;
                std::vector<double> colResults;
                std::vector<double> deltaResults;
                int rowIndex;
                int colIndex;
                int deltaIndex;
                int dataToGenerateNumber;
                std::vector<double> getRandomResults(std::vector<double>& results);
                double calcTotalEnergy();
                double getTotalEnergy();
                void calcInitialTotalEnergy();
                double calcAvgNumberOfSimilarNeighbours( int neighboursDistance, double limit );

                void setInitialData( double *data, int size );
                void setEnergyCalculator( EnergyCalculator *energyCalculator );
                void setRandomNumberGenerator( RandomNumberGenerator *randomNumberGenerator );
                void setMonterCarlo( MonteCarlo *mc );

                void setMaxChange( double maxChange );
                void setDataToChangeInSingleStep( int dataToChange );

                void singleStep();
                int getMaxNeighbours();

                Simulation();
};

#endif

