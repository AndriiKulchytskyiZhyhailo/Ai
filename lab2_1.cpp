#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <cstdlib>
#include <ctime>
#include <conio.h>

using namespace std;

double funkcja_test(double* koef, int n, double x);
double inth_funkcja_test(double* koef, int n, double a, double b);
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
    file_1.open("funkcja_test.txt");      // table funkcja test
    file_2.open("point__100_1.txt");
    file_3.open("point__100_0.txt");
    file_4.open("point_1000_1.txt");
    file_5.open("point_1000_0.txt");
    file_6.open("results.txt");
    double a = 0;
    double b = 1;
    int n_0 = 100;
    for (int i = 0; i <= n_0; i++)
    {
        double x = a + (b - a) * double(i) / double(n_0);
        file_1 << setw(12) << x << "   " << setw(12) << funkcja_test(koef_mn, n, x) << endl;
    }
    file_1.close();


    srand(131);
    //    srand(time(NULL));
    int N_point = 100;
    double f_max = 2.;
    double S_priam = (b - a) * f_max;
    double point[2];

    double inth_test = inth_funkcja_test(koef_mn, n, a, b);
    cout << "inth_funkcja_test = " << setw(10) << inth_test << endl << endl;
    file_6 << "inth_funkcja_test = " << setw(10) << inth_test << endl << endl;
    cout << " N_point      in_mn_cr          poch_abs      poch_vidn (%)" << endl << endl;
    file_6 << " N_point    in_mn_cr          poch_abs      poch_vidn (%)" << endl << endl;
    for (int k = 1; k <= 4; k++)
    {
        int n_point = 0;
        for (int i = 1; i <= N_point; i++)
        {
            point_vyp(a, b, f_max, point);

            if (point[1] <= funkcja_test(koef_mn, n, point[0]))
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
        double int_monte_carlo = S_priam * double(n_point) / double(N_point);
        double poch_abs = abs(inth_test - int_monte_carlo);
        double poch_vidn = fabs(inth_test - int_monte_carlo) * 100 / inth_test;
        cout << setw(8) << N_point << "    " << setw(10) << int_monte_carlo << "       ";
        cout << setw(12) << poch_abs << "    " << setw(12) << poch_vidn << endl;
        file_6 << setw(8) << N_point << "    " << setw(10) << int_monte_carlo << "       ";
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

double funkcja_test(double* koef, int n, double x)
{
    double mnohochlen = koef[0];
    for (int i = 1; i <= n; i++)
    {
        mnohochlen += koef[i] * x;
        x *= x;
    }
    return mnohochlen;
}

double inth_funkcja_test(double* koef, int n, double a, double b)
{
    double inth_mnohochlen = 0;
    for (int i = 0; i <= n; i++)
    {
        inth_mnohochlen += koef[i] * (b - a) / double(i + 1);
        a *= a;
        b *= b;
    }
    return inth_mnohochlen;
}