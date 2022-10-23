#ifndef SIMULATION_H

#define SIMULATION_H

#include<iostream>
#include"MyMPI.h"
#include"EnergyCalculator.h"
#include"RandomNumberGenerator.h"
#include"MonteCarlo.h"

using namespace std;

class Simulation {
private:
   
    MyMPI *mmpi;
    MonteCarlo *mc;
    RandomNumberGenerator *rng;
    EnergyCalculator *energyCalculator;
    int size;
    int reducedSize;
    int dataToChange;
    double maxChange;
    double *delta;
    int *rows;
    int *cols;
    double *data;
    double Etot;

    void changeData();
    void changeDataUndo();
    void generateDataChange();
    int calculateNumberOfRows(int sizeOfTable, int amountOfProcesses, int processRank);
    int calculateStartRow(int sizeOfTable, int amountOfProcesses, int processRank);

public:
    void init();
    double calcTotalEnergy(int rank, int numberOfProcesses,int startRow, int rowNumberEnd, int numberOfRows);
    double getTotalEnergy();
    void calcInitialTotalEnergy();

    void setInitialData( double *data, int size );
    void setEnergyCalculator( EnergyCalculator *energyCalculator );
    void setRandomNumberGenerator( RandomNumberGenerator *randomNumberGenerator );
    void setMonterCarlo( MonteCarlo *mc );

    void setMaxChange( double maxChange );
    void setDataToChangeInSingleStep( int dataToChange );

    void singleStep();

    Simulation( MyMPI *_mmpi );
};

#endif