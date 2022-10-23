#include <iostream>
#include <iomanip>


using namespace std;

void funkcja(double *table, int size){
    for(int i = 0; i < size; i++){
     cout<<i<<"   "<<table[i]<<endl;
    }
}


bool funkcja2(){
    return false;
}



int main(int ac, char **av){




    // for(int i = 0; i < size*size; i++){
    //     table[i] = i;
    // }

    // for(int i = 0; i < size*size; i++){
    //  cout<<i<<"   "<<table[i]<<endl;
    // }

    // int iloscProcesow = 3;

    // int iloscWierszy = (size)/3;

    // cout<<"Ilosc podstawowych wierszy: "<<iloscWierszy<<endl;
    // int iloscDodatkowychWierszy=size%iloscProcesow;
    // cout<<"Ilosc dodatkowych wierszy: "<<iloscDodatkowychWierszy<<endl;

    //proces 0 
    // int iloscWierszyDlaProcesu0 = iloscWierszy + 1;
    // int iloscPol = size * iloscWierszyDlaProcesu0;
    // cout<<"Ilosc Pol dla procesu 0: "<<iloscPol<<endl;

    // cout<<*(table+iloscPol)<<endl;
    // cout<<table<<endl;

    // cout<<"To teraz funkcja"<<endl;
    // funkcja(&table[size], size);


    bool flag = true;
    int wartosc = funkcja2();
    cout<<"Wartosc inta to: "<<wartosc<<endl;

    if(wartosc){
        cout<<"Dziala dla inta"<<endl;
    }else{
        cout<<"Nie dziala dla inta"<<endl;
    }



    // cout<<"Wielkosc tablicy: "<<sizeof(table)<<endl;



// int* wskaznik = &liczba;

// cout<<"Wartosc na ktora wskazuje wskaznik: "<<*wskaznik<<endl;
// cout<<"Adres na ktory wskazuje wskaznik: "<<wskaznik<<endl;
// cout<<"Adres wskaznika: "<<&wskaznik<<endl;




}




// Funkcja gdy sie udalo wyslac wskaznik
// void Simulation::init()
// {
//     int rank;
//     mmpi->MPI_Comm_rank(MPI_COMM_WORLD, &rank);

//     switch(rank){
//         case 0: {
//             char str[] = "napis do wyslania";
//             int liczba = 1000;
//             int* wskaznik = &liczba;
//             cout<<"wypisuje adres wskaznika: "<<&wskaznik<<endl;
//             cout<<"wypisuje adres wskaznika: "<<*wskaznik<<endl;
//             cout<<"wypisuje adres wskaznika: "<<wskaznik<<endl;


//         mmpi->MPI_Send(wskaznik, sizeof(int), MPI_INT, 1, 100, MPI_COMM_WORLD);
//         cout<<"Tu proces 0- dane wyslane"<<endl;
//         break;
//         }
//         case 1: {

//         cout<<"Tu proces 1 przed odebraniem danych"<<endl;
//         double* data2;
//         int sizeBuf = 512;
//         char *buf = new char[sizeBuf];
//         int liczba2;
//         int* wsk;
//         wsk = (int*)malloc(sizeof(int));


//         int wielkosc = sizeof(wsk);

//         MPI_Status status;
//         MPI_Recv(wsk, sizeof(int), MPI_INT, 0, 100, MPI_COMM_WORLD,&status);
//         // cout<<"tu proces 1- otrzymalem dane: "<<liczba2<<endl;
//         // cout<<"tu proces 1- otrzymalem dane: "<<&liczba2<<endl;
        
//         cout<<"tu proces 1- otrzymalem dane: "<<*wsk<<endl;
//         break;
//         }
//     }