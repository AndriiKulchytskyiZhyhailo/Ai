#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <cstdlib>
#include <ctime>
#include <conio.h>

using namespace std;

double funkcja_osnovna(double x);
void point_vyp(double a, double b, double h, double* point);

int main()
{
    int n = 2;
    double koef_mn[3] = { 1,4,-4 };

    ofstream file_1;
    ofstream file_2;
    ofstream file_3;
    ofstream file_4;
    ofstream file_5;
    ofstream file_6;
    file_1.open("funkcja_osnovna.txt");      // table funkcja test
    file_2.open("point__100_1_osn.txt");
    file_3.open("point__100_0_osn.txt");
    file_4.open("point_1000_1_osn.txt");
    file_5.open("point_1000_0_osn.txt");
    file_6.open("results_osn.txt");
    double a = 1;
    double b = 2;
    int n_0 = 100;
    for (int i = 0; i <= n_0; i++)
    {
        double x = a + (b - a) * double(i) / double(n_0);
        file_1 << setw(12) << x << "   " << setw(12) << funkcja_osnovna(x) << endl;
    }
    file_1.close();

    srand(131);
    //    srand(time(NULL));
    int N_point = 100;
    double f_max = exp(4.);
    double S_priam = (b - a) * f_max;
    double point[2];
    const int k_max = 6;
    double int_monte_carlo[k_max];

    for (int k = 1; k <= k_max; k++)
    {
        int n_point = 0;
        for (int i = 1; i <= N_point; i++)
        {
            point_vyp(a, b, f_max, point);

            if (point[1] <= funkcja_osnovna(point[0]))
            {
                n_point++;
                if (k == 1)
                    file_2 << setw(10) << point[0] << setw(10) << point[1] << endl;
                if (k == 2)
                    file_4 << setw(10) << point[0] << setw(10) << point[1] << endl;
            }
            else
            {
                if (k == 1)
                    file_3 << setw(10) << point[0] << setw(10) << point[1] << endl;
                if (k == 2)
                    file_5 << setw(10) << point[0] << setw(10) << point[1] << endl;
            }
        }
        int_monte_carlo[k - 1] = S_priam * double(n_point) / double(N_point);
        //cout << setw(8) << N_point << "    " << setw(10) << int_monte_carlo[k - 1] << endl;
        if (k == k_max)
        {
           // cout << endl << endl;
           // cout << "N_point = " << setw(8) << N_point << "   ";
            //cout << "inth_monte_carlo = " << setw(10) << int_monte_carlo[k - 1] << endl << endl;
            file_6 << "N_point = " << setw(8) << N_point << "   ";
            file_6 << "inth_monte_carlo = " << setw(10) << int_monte_carlo[k - 1] << endl << endl;
        }
        N_point *= 10;
    }
    cout << " N_point      in_mn_cr          poch_abs      poch_vidn (%)" << endl << endl;
    file_6 << " N_point    in_mn_cr          poch_abs      poch_vidn (%)" << endl << endl;

    N_point = 100;
    for (int k = 0; k < k_max-1 ; k++)
    {
        double poch_abs = abs(int_monte_carlo[k_max - 1] - int_monte_carlo[k]);
        double poch_vidn = fabs(int_monte_carlo[k_max - 1] - int_monte_carlo[k]) * 100 / int_monte_carlo[3];
        cout << setw(8) << N_point << "    " << setw(10) << int_monte_carlo[k] << "       ";
        cout << setw(12) << poch_abs << "    " << setw(12) << poch_vidn << endl;
        file_6 << setw(8) << N_point << "    " << setw(10) << int_monte_carlo[k] << "       ";
        file_6 << setw(12) << poch_abs << "    " << setw(12) << poch_vidn << endl;
        N_point *= 10;
    }
    
    return 0;
}

void point_vyp(double a, double b, double h, double* point)
{
    point[0] = a + (b - a) * double(rand()) / double(RAND_MAX);
    point[1] = h * double(rand()) / double(RAND_MAX);
}

double funkcja_osnovna(double x)
{
    return exp(x * x);
}
