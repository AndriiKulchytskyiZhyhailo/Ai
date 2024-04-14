#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <cstdlib>
#include <ctime>
#include <conio.h>

using namespace std;

void generator_test(double* point, double* center, int N, int n, double* center_test, int* n_point);
void oblicz_u(double* point, double* center, int N, int n, double* u, double par);
void modyfikacja(double* point, double* center, int N, int n, double* u, double par);
double miara_u(double* u_0, double* u_1, int N, int n);
void copy(double* u_0, double* u_1, int N, int n);
double suma_u(double* u, int N, int n);
double miara(double* x, double* c);
void vporiadkuvannia_vektor(double* x, int n, int* pozycia);

int main()
{
    const int N = 2000;  // number of points
    const int n = 5;     // number of centers 
    const int max_iter = 1000;
    double par = 2;
    double eps = 0.001;

    double point[2][N];
    double center[2][n];
    double u_0[n][N];
    double u_1[n][N];

    double center_test[3][n] = { {0.6,0.7,0.9,0.2,0.2},
                                  {0.3,0.7,0.1,0.8,0.2},
                                  {0.1,0.3,0.1,0.2,0.15} };
    int nm_points[n] = { 400,500,300,300,500 };
    generator_test(*point, *center, N, n, *center_test, nm_points);

    ofstream file_0_points;
    file_0_points.open("points_test_c.txt");
    for (int i = 0; i < N; i++)
        file_0_points << setw(12) << point[0][i] << "  " << setw(12) << point[1][i] << endl;
    file_0_points.close();

    ofstream file_0_centers;
    file_0_centers.open("centers_test_c.txt");
    for (int i = 0; i < n; i++)
    {
        file_0_centers << setw(10) << center_test[0][i] << "  " << setw(10) << center_test[1][i] << "  ";
        file_0_centers << setw(10) << center_test[2][i] << "  " << setw(10) << nm_points[i] << endl;
    }
    file_0_centers.close();

    ofstream file_1_centers;
    file_1_centers.open("centers_iteracie_c.txt");

    oblicz_u(*point, *center, N, n, *u_0, par);

    /////////////////////////////////////////////////////////////////////////////////////////////////////
    /////             iteracijnyj proces
    /////////////////////////////////////////////////////////////////////////////////////////////////////

    for (int iter = 0; iter < max_iter; iter++)
    {
        for (int j = 0; j < n; j++)
            file_1_centers << setw(12) << center[0][j] << "  " << setw(12) << center[1][j] << endl;
        file_1_centers << endl;
        modyfikacja(*point, *center, N, n, *u_0, par);
        oblicz_u(*point, *center, N, n, *u_1, par);
        double norma = miara_u(*u_0, *u_1, N, n);
        cout << "nr_iter = " << setw(4) << iter << "  norma= " << setw(6) << norma << endl;
        if (norma < eps) break;
        copy(*u_0, *u_1, N, n);
    }
    file_1_centers.close();

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////      vypysuvannia rezultativ
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    ofstream file_1_points_klasters;
    file_1_points_klasters.open("points_klasters_c1.txt");

    double* stovpchyk_u = new double[n];
    int* i_stovpchyk_u = new int[n];
    int* point_nr_center = new int[N];
    double* point_p_center = new double[N];
    double limit = 0.5;
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < n; j++)
            stovpchyk_u[j] = u_1[j][i];
        vporiadkuvannia_vektor(stovpchyk_u, n, i_stovpchyk_u);
        if (stovpchyk_u[0] > limit)
        {
            point_nr_center[i] = i_stovpchyk_u[0] + 1;
            point_p_center[i] = stovpchyk_u[0];
        }
        else
        {
            point_nr_center[i] = (i_stovpchyk_u[0] + 1) * 10 + i_stovpchyk_u[1] + 1;
            point_p_center[i] = stovpchyk_u[0] + stovpchyk_u[1];
        }
        if (point_p_center[i] < limit)
        {
            cout << "point nr " << setw(4) << i << setw(10) << point_nr_center[i] << endl;
        }
    }
    for (int k = 0; k < n; k++)
        for (int i = 0; i < N; i++)
            if (point_nr_center[i] - 1 == k)
            {
                file_1_points_klasters << setw(12) << point[0][i];
                file_1_points_klasters << setw(12) << point[1][i];
                file_1_points_klasters << setw(12) << point_nr_center[i];
                file_1_points_klasters << setw(12) << point_p_center[i] << endl;
            }
    for (int i = 0; i < N; i++)
        if (point_nr_center[i] > 9)
        {
            file_1_points_klasters << setw(12) << point[0][i];
            file_1_points_klasters << setw(12) << point[1][i];
            file_1_points_klasters << setw(12) << point_nr_center[i];
            file_1_points_klasters << setw(12) << point_p_center[i] << endl;
        }

    file_1_points_klasters.close();

    ofstream file_2_points_klasters;
    file_2_points_klasters.open("points_klasters_c2.txt");

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < n; j++)
            stovpchyk_u[j] = u_1[j][i];
        vporiadkuvannia_vektor(stovpchyk_u, n, i_stovpchyk_u);
        point_nr_center[i] = int(stovpchyk_u[0] / 0.1);
        point_p_center[i] = stovpchyk_u[0];
    }
    for (int k = 9; k >= 0; k--)
        for (int i = 0; i < N; i++)
            if (point_nr_center[i] - 1 == k)
            {
                file_2_points_klasters << setw(12) << point[0][i];
                file_2_points_klasters << setw(12) << point[1][i];
                file_2_points_klasters << setw(12) << point_nr_center[i];
                file_2_points_klasters << setw(12) << point_p_center[i] << endl;
            }
    file_2_points_klasters.close();

    return 0;
}

void oblicz_u(double* point, double* center, int N, int n, double* u, double par)
{
    double* temp_array = new double[n];
    for (int i = 0; i < N; i++)
    {
        double temp = 0;
        for (int j = 0; j < n; j++)
        {
            double x[2];
            double c[2];
            x[0] = point[0 * N + i];
            x[1] = point[1 * N + i];
            c[0] = center[0 * n + j];
            c[1] = center[1 * n + j];
            temp_array[j] = exp(log(miara(x, c)) * 2. / (1 - par));
            temp += temp_array[j];
        }
        for (int j = 0; j < n; j++)
            u[j * N + i] = temp_array[j] / temp;
    }
}

void modyfikacja(double* point, double* center, int N, int n, double* u, double par)
{
    for (int j = 0; j < n; j++)
    {
        double temp_x = 0;
        double temp_y = 0;
        double temp_0 = 0;
        for (int i = 0; i < N; i++)
        {
            temp_x += pow(u[j * N + i], par) * point[0 * N + i];
            temp_y += pow(u[j * N + i], par) * point[1 * N + i];
            temp_0 += pow(u[j * N + i], par);
        }
        center[0 * n + j] = temp_x / temp_0;
        center[1 * n + j] = temp_y / temp_0;
    }
}

double miara(double* x, double* c)
{
    return sqrt(pow(x[0] - c[0], 2) + pow(x[1] - c[1], 2));
}

double miara_u(double* u_0, double* u_1, int N, int n)
{
    double miara_m = 0;
    for (int i = 0; i < N * n; i++)
        miara_m += pow(u_0[i] - u_1[i], 2);
    return sqrt(miara_m);
}

double suma_u(double* u, int N, int n)
{
    double suma = 0;
    for (int i = 0; i < N * n; i++)
        suma += u[i];
    return suma;
}

void copy(double* u_0, double* u_1, int N, int n)
{
    for (int i = 0; i < N * n; i++)
        u_0[i] = u_1[i];
}

void generator_test(double* point, double* center, int N, int n, double* center_test, int* n_point)
{
    srand(131);
    int nr = 0;
    int* nr_k = new int[n];
    for (int k = 0; k < n; k++)
        nr_k[k] = 0;
    do
    {
        double x = double(rand()) / double(RAND_MAX);
        double y = double(rand()) / double(RAND_MAX);
        for (int k = 0; k < n; k++)
        {
            double temp = pow(x - center_test[0 * n + k], 2) + pow(y - center_test[1 * n + k], 2);
            if ((temp < pow(center_test[2 * n + k], 2)) && (nr_k[k] < n_point[k]))
            {
                point[0 * N + nr] = x;
                point[1 * N + nr] = y;
                nr++;
                nr_k[k]++;
                break;
            }
        }
    } while (nr < N);

    //    for (int i=0;i<n;i++)
    //    for (int k=0;k<2;k++)
    //        center[k*n+i] = double(rand())/double(RAND_MAX);

    center[0] = 0.50;
    center[5] = 0.50;
    center[1] = 0.25;
    center[6] = 0.25;
    center[2] = 0.25;
    center[7] = 0.75;
    center[3] = 0.75;
    center[8] = 0.25;
    center[4] = 0.75;
    center[9] = 0.75;
}

void vporiadkuvannia_vektor(double* x, int n, int* pozycia)
{
    for (int i = 0; i < n; i++)
        pozycia[i] = i;
    for (int i = n - 1; i > 0; i--)
        for (int j = 0; j < i; j++)
            if (x[j] < x[j + 1])
            {
                double temp = x[j];
                x[j] = x[j + 1];
                x[j + 1] = temp;
                int i_temp = pozycia[j];
                pozycia[j] = pozycia[j + 1];
                pozycia[j + 1] = i_temp;
            }
}