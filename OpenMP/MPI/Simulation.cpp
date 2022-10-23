#include "Simulation.h"
#include "Helper.h"
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <iomanip>

using namespace std;

Simulation::Simulation(MyMPI *_mmpi)
{
    mmpi = _mmpi;
}

void Simulation::setRandomNumberGenerator(RandomNumberGenerator *randomNumberGenerator)
{
    this->rng = randomNumberGenerator;
}

void Simulation::setEnergyCalculator(EnergyCalculator *energyCalculator)
{
    this->energyCalculator = energyCalculator;
}

void Simulation::setMonterCarlo(MonteCarlo *mc)
{
    this->mc = mc;
}

void Simulation::setMaxChange(double maxChange)
{
    this->maxChange = maxChange;
}

void Simulation::setInitialData(double *data, int size)
{
    this->data = data;
    this->size = size;
    reducedSize = size - 4;
}

void Simulation::init()
{

    int rank;
    int numberOfProcesses;
    mmpi->MPI_Comm_size(MPI_COMM_WORLD, &numberOfProcesses);
    mmpi->MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (numberOfProcesses > 1)
    {
        int position = 0;

        switch (rank)
        {
        case 0:
        {
            mmpi->MPI_Bcast(&size, sizeof(int), MPI_INT, 0, MPI_COMM_WORLD);
            mmpi->MPI_Bcast(data, size * size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            mmpi->MPI_Bcast(&this->dataToChange, sizeof(int), MPI_INT, 0, MPI_COMM_WORLD);
        }
        break;
        default:
        {
            mmpi->MPI_Bcast(&size, sizeof(int), MPI_INT, 0, MPI_COMM_WORLD);
            data = new double[size * size];
            mmpi->MPI_Bcast(data, size * size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            mmpi->MPI_Bcast(&this->dataToChange, sizeof(int), MPI_INT, 0, MPI_COMM_WORLD);
            rows = new int[dataToChange];
            cols = new int[dataToChange];
            delta = new double[dataToChange];
        }
            break;
        }
    }

}

void Simulation::calcInitialTotalEnergy()
{

    int rank;
    int numberOfProcesses;
    int numberOfRows;
    int startRow;
    int rowNumberEnd;
    double momentaryEtot = 0.0;

    mmpi->MPI_Comm_size(MPI_COMM_WORLD, &numberOfProcesses);
    mmpi->MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    numberOfRows = calculateNumberOfRows(size, numberOfProcesses, rank);
    startRow = calculateStartRow(size, numberOfProcesses, rank);
    rowNumberEnd = startRow + numberOfRows;
    momentaryEtot = this->calcTotalEnergy(rank, numberOfProcesses, startRow, rowNumberEnd, numberOfRows);
    mmpi->MPI_Reduce(&momentaryEtot, &Etot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
}

double Simulation::calcTotalEnergy(int rank, int numberOfProcesses,int startRow, int rowNumberEnd, int numberOfRows)
{
    double Etot = 0.0;

    if (rank == 0)
        startRow = 2;

    if (rank == numberOfProcesses - 1)
        rowNumberEnd -= 2;

    for (startRow; startRow < rowNumberEnd; startRow++)
        for (int col = 2; col < size - 2; col++)
            Etot += energyCalculator->calc(data, size, startRow, col);

    return Etot * 0.5;
}


double Simulation::getTotalEnergy()
{
    return Etot;
}

void Simulation::setDataToChangeInSingleStep(int dataToChange)
{
    this->dataToChange = dataToChange;
    rows = new int[dataToChange];
    cols = new int[dataToChange];
    delta = new double[dataToChange];
}

void Simulation::generateDataChange()
{
    for (int i = 0; i < dataToChange; i++)
    {
        rows[i] = 2 + rng->getInt(reducedSize);
        cols[i] = 2 + rng->getInt(reducedSize);
        delta[i] = maxChange * (1.0 - 2.0 * rng->getDouble());
    }
}

void Simulation::changeData()
{
    for (int i = 0; i < dataToChange; i++)
    {
        Helper::updateValue(data, size, rows[i], cols[i], delta[i]);
    }
}

void Simulation::changeDataUndo()
{
    for (int i = 0; i < dataToChange; i++)
    {
        Helper::updateValue(data, size, rows[i], cols[i], -delta[i]);
    }
}

void Simulation::singleStep()
{

    int rank;
    int numberOfProcesses;
    int flag;
    int numberOfRows;
    int startRow;
    int rowNumberEnd;
    double momentaryEtot = 0.0;
    MPI_Status status;

    mmpi->MPI_Comm_size(MPI_COMM_WORLD, &numberOfProcesses);
    mmpi->MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    numberOfRows = calculateNumberOfRows(size, numberOfProcesses, rank);
    startRow = calculateStartRow(size, numberOfProcesses, rank);
    rowNumberEnd = startRow + numberOfRows;

    switch (rank)
    {
    case 0:
    {
        generateDataChange();
        for (int i = 1; i < numberOfProcesses; i++)
        {
            mmpi->MPI_Send(rows, dataToChange, MPI_INT, i, i, MPI_COMM_WORLD);
            mmpi->MPI_Send(cols, dataToChange, MPI_INT, i, i, MPI_COMM_WORLD);
            mmpi->MPI_Send(delta, dataToChange, MPI_DOUBLE, i, i, MPI_COMM_WORLD);
        }
    }
    break;
    default:
    {
        
        mmpi->MPI_Recv(rows, dataToChange, MPI_INT, 0, rank, MPI_COMM_WORLD, &status);
        mmpi->MPI_Recv(cols, dataToChange, MPI_INT, 0, rank, MPI_COMM_WORLD, &status);
        mmpi->MPI_Recv(delta, dataToChange, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD, &status);

    }
    break;
    }

    changeData(); // zmiana danych (stanu)
    momentaryEtot = calcTotalEnergy(rank, numberOfProcesses, startRow, rowNumberEnd, numberOfRows);
    double newEtot = 0.0;
    mmpi->MPI_Reduce(&momentaryEtot, &newEtot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    switch (rank)
    {
    case 0:
    {
        flag = mc->accept(Etot, newEtot);
        for (int i = 1; i < numberOfProcesses; i++)
        {
            mmpi->MPI_Send(&flag, 1, MPI_INT, i, i, MPI_COMM_WORLD);

        }
    }
    break;
    default:
    {
        mmpi->MPI_Recv(&flag, 1, MPI_INT, 0, rank, MPI_COMM_WORLD, &status);
    }

    break;
    }


    if (flag)
    {
        if(rank == 0){
            cout << "Accepted Eold " << Etot << " newE " << newEtot << endl;
        }
        Etot = newEtot;
    }else
    {
        changeDataUndo();
        if(rank == 0){
            cout << "Not accepted Eold " << Etot << " newE " << newEtot << endl;
        }
    }
}

int Simulation::calculateNumberOfRows(int sizeOfTable, int numberOfProcesses, int processRank)
{
    int extraAmountOfRows = sizeOfTable % numberOfProcesses;
    int basicAmountOfRows = sizeOfTable / numberOfProcesses;
    int difference = numberOfProcesses - extraAmountOfRows;

    if (extraAmountOfRows == 0)
    {
        return basicAmountOfRows;
    }else
    {
        if (processRank >= difference)
        {
            return basicAmountOfRows + 1;
        }
        else
        {
            return basicAmountOfRows;
        }
    }
}

int Simulation::calculateStartRow(int sizeOfTable, int numberOfProcesses, int processRank)
{
    int startRow = 0;

    for (int i = 0; i < processRank; i++)
    {
        startRow += calculateNumberOfRows(sizeOfTable, numberOfProcesses, i);
    }
    return startRow;
}