#include "Simulation.h"
#include "Helper.h"
#include <iostream>
#include <math.h>
#include <iomanip>
#include <stdio.h>
#include <omp.h>

using namespace std;

// Simulation::Simulation()
// {
// }

Simulation::Simulation()
{
#pragma omp parallel
    {
        if (omp_get_thread_num() == 0)
        {
            num_threads = omp_get_num_threads();
            printf("The number of threads: %d\n", num_threads);
            // cout << "Number of threads in constructor: " << num_threads << endl;
            this->structures = new drand48_data[num_threads];
        }
    }

    #pragma omp parallel for
    for ( int i = 0; i < omp_get_num_threads(); i++){
        srand48_r(omp_get_thread_num(), &structures[omp_get_thread_num()]);
        // double res;
        // drand48_r( &structures[i], &res);
        // printf( "Example: %f\n", res );
    }
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

void Simulation::calcInitialTotalEnergy()
{
   Etot = this->calcTotalEnergy();
}

// do zrĂłwnoleglenia
double Simulation::calcTotalEnergy()
{
   double Etot = 0.0;
   double r_Etot = 0.0;
   // int id = omp_get_thread_num();
   // int nthreads = omp_get_num_threads();
   //
   // printf("This thread: %d, all threads: %d\n", id, nthreads);

#pragma omp parallel for reduction(+: Etot)
   for (int row = 2; row < size - 2; row++)
      for (int col = 2; col < size - 2; col++)
         Etot += energyCalculator->calc(data, size, row, col);
   // # pragma omp parallel reduction( +:Etot )
   // {
   //     int lines_per_process = floor( size / omp_get_num_threads());
   //     // int rest = size - omp_get_num_threads() * lines_per_process;
   //     int start = lines_per_process * omp_get_thread_num();
   //     int end = start + lines_per_process;
   //     if (omp_get_thread_num() == 0)
   //     {
   //         start = 2;
   //         lines_per_process -= 2;
   //     }
   //     else if (omp_get_thread_num() == omp_get_num_threads()-1)
   //     {
   //         end = size - 2;
   //         lines_per_process -= 2;
   //     }
   //    for (int row = start; row < end; row++)
   //       for (int col = 2; col < size - 2; col++)
   //          Etot += energyCalculator->calc(data, size, row, col);
   // }

   return Etot * 0.5;
}

int Simulation::similarNeighbours(int col, int row, int delta, double limit)
{
   double middle = Helper::getValue(data, size, row, col);
   int neighbours = 0;
   for (int rowd = -delta; rowd <= delta; rowd++)
      for (int cold = -delta; cold <= delta; cold++)
      {
         if (cos(Helper::getValue(data, size, row + rowd, col + cold) - middle) > limit)
            neighbours++;
      }

   return neighbours - 1;
}

// do zrĂłwnoleglenia
double Simulation::calcAvgNumberOfSimilarNeighbours(int neighboursDistance, double limit)
{
   int sum = 0;
   int neighbours = 0;
   maxNeighbours = 0;
   int neighboursTmp;
   int localMax = 0;
#pragma omp parallel for reduction(+ \
                                   : neighbours, sum) firstprivate(neighboursTmp, localMax)
   for (int row = neighboursDistance; row < size - neighboursDistance; row++)
   {
      for (int col = neighboursDistance; col < size - neighboursDistance; col++)
      {
         neighboursTmp = similarNeighbours(col, row, neighboursDistance, limit);
         sum += neighboursTmp;
         if (neighboursTmp > localMax)
         {
            localMax = neighboursTmp;
         }
         neighbours++;
      }
#pragma omp critical
      if (localMax > maxNeighbours)
      {
         maxNeighbours = localMax;
      }
   }
   return (double)sum / (double)(neighbours * (neighboursDistance + 1) * 4 * neighboursDistance);
}

double Simulation::getTotalEnergy()
{
   return Etot;
}

int Simulation::getMaxNeighbours()
{
   return maxNeighbours;
}

void Simulation::setDataToChangeInSingleStep(int dataToChange)
{
   this->dataToChange = dataToChange;
   rows = new int[dataToChange];
   cols = new int[dataToChange];
   delta = new double[dataToChange];
}

// do zrĂłwnoleglenia - konieczna wymiara generatora liczb losowych na drand48_r
void Simulation::generateDataChange()
{
   double d_res = 0.0;
#pragma omp parallel for firstprivate(d_res)
   for (int i = 0; i < dataToChange; i++)
   {
      // rows[i] = 2 + rng->getInt(reducedSize);
      // cols[i] = 2 + rng->getInt(reducedSize);
      // delta[i] = maxChange * (1.0 - 2.0 * rng->getDouble());
      drand48_r(&structures[omp_get_thread_num()], &d_res);
      rows[i] = 2 + (int)(d_res * reducedSize);
      drand48_r(&structures[omp_get_thread_num()], &d_res);
      cols[i] = 2 + (int)(d_res * reducedSize);
      drand48_r(&structures[omp_get_thread_num()], &d_res);
      delta[i] = maxChange * (1.0 - 2.0 * d_res);
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
   generateDataChange(); // wygenerowanie danych potrzebnych do zmiany stanu
   changeData();         // zmiana danych (stanu)

   // calcTotalEnergy
   double newEtot = calcTotalEnergy(); // wyliczenie nowej energii caĹkowitej

   // decyzja modulu MonteCarlo o akceptacji zmiany
   if (mc->accept(Etot, newEtot))
   {
      cout << "Accepted Eold " << Etot << " newE " << newEtot << endl;
      Etot = newEtot;
      // zaakceptowano zmiane -> nowa wartosc energii calkowitej
   }
   else
   {
      changeDataUndo();
      cout << "Not accepted Eold " << Etot << " newE " << newEtot << endl;
      // zmiany nie zaakceptowano -> przywracany stary stan, energia bez zmiany
   }
}
