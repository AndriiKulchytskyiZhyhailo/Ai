#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <cstdlib>
#include <ctime>
#include <conio.h>
#include <chrono> 

using namespace std;

const int N = 8;
short int* solution;
short int* columna;
short int* diagonal_1;
short int* diagonal_2;
short int* start_riadok;

ofstream file;

void inp_outp_queen(short int i, short int j, short int k);
void print_solution();
void fprint_solution();

int main()
{

    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;



    file.open("FERZI_SZACHOVA_DOSHKA.txt");      // table funkcja test

    solution = new short int[N];                // rozviazok 
    columna = new short int[N];
    diagonal_1 = new short int[2 * N - 1];
    diagonal_2 = new short int[2 * N - 1];
    start_riadok = new short int[N];
    

    for (short int i = 0; i < N; i++)
    {
        columna[i] = 1;                                
        start_riadok[i] = 0;              
    }                                     
    for (short int i = 0; i < 2 * N - 1; i++)
    {
        diagonal_1[i] = 1;                
        diagonal_2[i] = 1;                
    }
    
    //srand(time(NULL));
    srand(131);                              // % ostacza vid dilennia: rezultat vid 0 do 7
    short int pozycia_start = rand() % 8;               
    cout << "pozycia startowa = " << " 1  " << char(pozycia_start + 65) << endl; // znak A maje kod 65; znak B - 66
    file << "pozycia startowa = " << " 1  " << char(pozycia_start + 65) << endl; // i tak dali
    // 
    solution[0] = pozycia_start;
    inp_outp_queen(0, pozycia_start, 0);   
    // pewni kolonky i diahonali (vylajemo tretij arhument 0)
    auto t1 = high_resolution_clock::now();
    for (short int i = 1; i < N; i++)          
    {

        //        cout << "i= " << i+1 ;
        short int j;
        for (j = start_riadok[i]; j < N; j++)     
            if (columna[j] && diagonal_1[i + j] && diagonal_2[N - 1 + i - j])  // 
                break;                      
        //        cout << " j= " << j+1 << endl;    // ne znajszly: j=8 (normalne zkinchennia cyklu, ne break
        if (j < N)
        {                            
            solution[i] = j;
            if (i + 1 == N) break;     
            inp_outp_queen(i, j, 0);   
            start_riadok[i] = j + 1;  
            //            cout << "i= " << setw(1) << i+1 << "    j= " << setw(1) << j+1 << endl;
        }
        else                         
        {
            start_riadok[i] = 0;           
            i--;                     
            inp_outp_queen(i, solution[i], 1);  
            i--;                     
        }
    }
    auto t2 = high_resolution_clock::now();
    auto ms_int = duration_cast<milliseconds>(t2 - t1);
    cout << endl << endl << "                         SOLUTION" << endl << endl << endl;
    file << endl << endl << "                         SOLUTION" << endl << endl << endl;
    
    print_solution();
    fprint_solution();
    cout <<"\n\ntime_of working = " << ms_int.count()<<endl;

   // cout << "Time taken: " << t.elapsed() << '\n';

    file.close();
    return 0;
}

void inp_outp_queen(short int i, short int j, short int k)
{
    columna[j] = k;
    diagonal_1[i + j] = k;
    diagonal_2[N - 1 + i - j] = k;
}

void print_solution()
{
    cout << "      ";
    for (short int i = 0; i < N; i++)
        cout << "  " << char(65 + i) << "   ";
    cout << endl;
    cout << "     _________________________________________________" << endl;
    for (short int j = N; j > 0; j--)
    {
        cout << "     |     |     |     |     |     |     |     |     |" << endl;
        cout << "  " << j << "  |";
        for (short int k = 0; k < N; k++)
        {
            if (k == solution[j - 1])
                cout << " FFF " << "|";
            else
                cout << "     " << "|";
        }
        cout << "  " << j << endl;
        cout << "     |     |     |     |     |     |     |     |     |" << endl;
        cout << "     _________________________________________________" << endl;
    }
    cout << endl << "      ";
    for (short int i = 0; i < N; i++)
        cout << "  " << char(65 + i) << "   ";
}

void fprint_solution()
{
    file << "      ";
    for (short int i = 0; i < N; i++)
        file << "  " << char(65 + i) << "   ";
    file << endl;
    file << "     _________________________________________________" << endl;
    for (short int j = N; j > 0; j--)
    {
        file << "     |     |     |     |     |     |     |     |     |" << endl;
        file << "  " << j << "  |";
        for (short int k = 0; k < N; k++)
        {
            if (k == solution[j - 1])
                file << " FFF " << "|";
            else
                file << "     " << "|";
        }
        file << "  " << j << endl;
        file << "     |     |     |     |     |     |     |     |     |" << endl;
        file << "     _________________________________________________" << endl;
    }
    file << endl << "      ";
    for (short int i = 0; i < N; i++)
        file << "  " << char(65 + i) << "   ";
}