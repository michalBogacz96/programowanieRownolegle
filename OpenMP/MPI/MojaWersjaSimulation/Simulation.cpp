#include "Simulation.h"
#include "Helper.h"
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <iomanip>

using namespace std;

Simulation::Simulation(MyMPI *_mmpi){
    mmpi = _mmpi;
}

void Simulation::setRandomNumberGenerator(RandomNumberGenerator *randomNumberGenerator){
    this->rng = randomNumberGenerator;
}

void Simulation::setEnergyCalculator(EnergyCalculator *energyCalculator){
    this->energyCalculator = energyCalculator;
}

void Simulation::setMonterCarlo(MonteCarlo *mc)
{
    this->mc = mc;
}

void Simulation::setMaxChange(double maxChange){
    this->maxChange = maxChange;
}

void Simulation::setInitialData(double *data, int size){
    this->data = data;
    this->size = size;
    reducedSize = size - 4;
}

void Simulation::init(){

    int rank;
    int amountOfProcesses;
    mmpi->MPI_Comm_size(MPI_COMM_WORLD, &amountOfProcesses);
    mmpi->MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (amountOfProcesses > 1){
        int position = 0;

        switch (rank){
        case 0:
        {
            mmpi->MPI_Bcast(&size, 1, MPI_INT, 0, MPI_COMM_WORLD);
            mmpi->MPI_Bcast(data, size * size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            mmpi->MPI_Bcast(&this->dataToChange, 1, MPI_INT, 0, MPI_COMM_WORLD);
        }
        break;
        default:
        {

            mmpi->MPI_Bcast(&size, 1, MPI_INT, 0, MPI_COMM_WORLD);
            data = new double[size * size];
            mmpi->MPI_Bcast(data, size * size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            mmpi->MPI_Bcast(&this->dataToChange, 1, MPI_INT, 0, MPI_COMM_WORLD);
            rows = new int[dataToChange];
            cols = new int[dataToChange];
            delta = new double[dataToChange];
        }
            break;
        }
    }

    this->startRow = calculateStartRow(size, amountOfProcesses, rank);
    this->numberOfRows = calculateNumberOfRows(size, amountOfProcesses, rank);
}

int Simulation::calculateNumberOfRows(int sizeOfTable, int amountOfProcesses, int processRank){
    int basicAmountOfRows = sizeOfTable / amountOfProcesses;
    int extraAmountOfRows = sizeOfTable % amountOfProcesses;
    int amountOfRows = -1;

    if (extraAmountOfRows == 0){
        amountOfRows = basicAmountOfRows;
    }
    else{
        if (processRank < extraAmountOfRows){
            amountOfRows = basicAmountOfRows + 1;
        }
        else{
            amountOfRows = basicAmountOfRows;
        }
    }
    return amountOfRows;
}

int Simulation::calculateStartRow(int sizeOfTable, int amountOfProcesses, int processRank){
    int startingRow = 0;

    for (int i = 0; i < processRank; i++){
        startingRow += calculateNumberOfRows(sizeOfTable, amountOfProcesses, i);
    }

        if (processRank == 0){
            startingRow = 2;
        }

    return startingRow;
}

void Simulation::calcInitialTotalEnergy(){

    int rank;
    int amountOfProcesses;
    mmpi->MPI_Comm_size(MPI_COMM_WORLD, &amountOfProcesses);
    mmpi->MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (amountOfProcesses > 1){
        double tmpEtot = 0.0;
        tmpEtot = this->calcTotalEnergy(rank, amountOfProcesses);
        mmpi->MPI_Reduce(&tmpEtot, &Etot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }
    else{
        Etot = this->calcTotalEnergy(rank, amountOfProcesses);
    }
}

double Simulation::calcTotalEnergy(int rank, int amountOfProcesses){
    double Etot = 0.0;
    int endRow = startRow + numberOfRows;

    if (rank == 0 || rank == amountOfProcesses - 1){
        endRow -= 2;
    }

    for (int i = startRow; i < endRow; i++){
        for (int col = 2; col < size - 2; col++){
            Etot += energyCalculator->calc(data, size, i, col);
        }
    }
    return Etot * 0.5;
}


double Simulation::getTotalEnergy(){
    return Etot;
}

void Simulation::setDataToChangeInSingleStep(int dataToChange){
    this->dataToChange = dataToChange;
    rows = new int[dataToChange];
    cols = new int[dataToChange];
    delta = new double[dataToChange];
}

void Simulation::generateDataChange(){
    for (int i = 0; i < dataToChange; i++){
        rows[i] = 2 + rng->getInt(reducedSize);
        cols[i] = 2 + rng->getInt(reducedSize);
        delta[i] = maxChange * (1.0 - 2.0 * rng->getDouble());
    }
}

void Simulation::changeData(){
    for (int i = 0; i < dataToChange; i++){
        Helper::updateValue(data, size, rows[i], cols[i], delta[i]);
    }
}

void Simulation::changeDataUndo(){
    for (int i = 0; i < dataToChange; i++){
        Helper::updateValue(data, size, rows[i], cols[i], -delta[i]);
    }
}

void Simulation::singleStep(){
    
    int rank;
    int amountOfProcesses;
    int flag;
    double newEtot = 0.0;
    mmpi->MPI_Comm_size(MPI_COMM_WORLD, &amountOfProcesses);
    mmpi->MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    switch (rank){
    case 0:
    {
        generateDataChange();
        mmpi->MPI_Bcast(rows, dataToChange, MPI_INT, 0, MPI_COMM_WORLD);
        mmpi->MPI_Bcast(cols, dataToChange, MPI_INT, 0, MPI_COMM_WORLD);
        mmpi->MPI_Bcast(delta, dataToChange, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
    break;
    default:
    {
        mmpi->MPI_Bcast(rows, dataToChange, MPI_INT, 0, MPI_COMM_WORLD);
        mmpi->MPI_Bcast(cols, dataToChange, MPI_INT, 0, MPI_COMM_WORLD);
        mmpi->MPI_Bcast(delta, dataToChange, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    break;
    }

    changeData();
    double tmpEtot = calcTotalEnergy(rank, amountOfProcesses);

    mmpi->MPI_Reduce(&tmpEtot, &newEtot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    switch (rank){
    case 0:
    {
        flag = mc->accept(Etot, newEtot);
        mmpi->MPI_Bcast(&flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }
    break;
    default:
    {
        mmpi->MPI_Bcast(&flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }

    break;
    }

    if (flag){
        if(rank == 0){
            cout << "Accepted Eold " << Etot << " newE " << newEtot << endl;
        }
        Etot = newEtot;
    }else{
        changeDataUndo();
        if(rank == 0){
            cout << "Not accepted Eold " << Etot << " newE " << newEtot << endl;
        }
    }
}