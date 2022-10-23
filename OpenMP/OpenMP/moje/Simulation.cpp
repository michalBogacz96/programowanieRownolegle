#include "Simulation.h"
#include "Helper.h"
#include <iostream>
#include <math.h>
#include <iomanip>
#include <stdio.h>
#include <omp.h>

using namespace std;
#define MY_CONST 64

// Simulation::Simulation()
// {
// }

Simulation::Simulation()
{

   int numberOfAllThreads = 0;
   #pragma omp parallel firstprivate(numberOfAllThreads)
   {
      numberOfAllThreads = omp_get_num_threads();
      if (omp_get_thread_num() == 0)
      {
         this->generator = new drand48_data[MY_CONST * numberOfAllThreads];
      }
   }

   int numberOfTheThread = 0;
   #pragma omp parallel for firstprivate(numberOfTheThread)
   for (int i = 0; i < omp_get_num_threads(); i++)
   {
      numberOfTheThread = omp_get_thread_num();
      srand48_r(numberOfTheThread, &generator[MY_CONST * numberOfTheThread]);
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

//// do zrĂłwnoleglenia
//double Simulation::calcTotalEnergy()
//{
//   double Etot = 0.0;
//   for (int row = 2; row < size - 2; row++)
//      for (int col = 2; col < size - 2; col++)
//         Etot += energyCalculator->calc(data, size, row, col);
//   return Etot * 0.5;
//}

// do zrĂłwnoleglenia
double Simulation::calcTotalEnergy(){
   double Etot = 0.0;

   #pragma omp parallel for reduction(+: Etot)
   for (int row = 2; row < size - 2; row++)
      for (int col = 2; col < size - 2; col++)
         Etot += energyCalculator->calc(data, size, row, col);
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
//double Simulation::calcAvgNumberOfSimilarNeighbours(int neighboursDistance, double limit)
//{
//   int sum = 0;
//   int neighbours = 0;
//   maxNeighbours = 0;
//   int neighboursTmp;
//   for (int row = neighboursDistance; row < size - neighboursDistance; row++)
//      for (int col = neighboursDistance; col < size - neighboursDistance; col++)
//      {
//         neighboursTmp = similarNeighbours(col, row, neighboursDistance, limit);
//         sum += neighboursTmp;
//         if (neighboursTmp > maxNeighbours)
//         {
//            maxNeighbours = neighboursTmp;
//         }
//         neighbours++;
//      }
//   return (double)sum / (double)(neighbours * (neighboursDistance + 1) * 4 * neighboursDistance);
//}

// do zrĂłwnoleglenia
double Simulation::calcAvgNumberOfSimilarNeighbours(int neighboursDistance, double limit)
{
   int maxNeighboursTmp = 0;
   int sum = 0;
   int neighbours = 0;
   maxNeighbours = 0;
   int neighboursTmp;
   

   #pragma omp parallel for reduction(+: sum, neighbours) firstprivate(neighboursTmp, maxNeighboursTmp)
   for (int row = neighboursDistance; row < size - neighboursDistance; row++){

      for (int col = neighboursDistance; col < size - neighboursDistance; col++){
         neighboursTmp = similarNeighbours(col, row, neighboursDistance, limit);
         sum += neighboursTmp;

         if (neighboursTmp > maxNeighboursTmp)
         {
            maxNeighboursTmp = neighboursTmp;
         }

         neighbours++;
      }

      #pragma omp critical
      if (maxNeighbours < maxNeighboursTmp)
      {
         maxNeighbours = maxNeighboursTmp;
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
// void Simulation::generateDataChange()
// {
//    for (int i = 0; i < dataToChange; i++)
//    {
//       rows[i] = 2 + rng->getInt(reducedSize);
//       cols[i] = 2 + rng->getInt(reducedSize);
//       delta[i] = maxChange * (1.0 - 2.0 * rng->getDouble());
//    }
// }

// do zrĂłwnoleglenia - konieczna wymiara generatora liczb losowych na drand48_r
void Simulation::generateDataChange()
{
   double generatedNumber = 0.0;
   int numberOfTheThread = 0;

   #pragma omp parallel for firstprivate(generatedNumber, numberOfTheThread)
   for (int i = 0; i < dataToChange; i++)
   {
      numberOfTheThread = omp_get_thread_num();
      drand48_r(&generator[MY_CONST * numberOfTheThread], &generatedNumber);
      rows[i] = 2 + (int)(reducedSize * generatedNumber);
      drand48_r(&generator[MY_CONST * numberOfTheThread], &generatedNumber);
      cols[i] = 2 + (int)(reducedSize * generatedNumber);
      drand48_r(&generator[MY_CONST * numberOfTheThread], &generatedNumber);
      delta[i] = maxChange * (1.0 - 2.0 * generatedNumber);
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
