#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <cstdlib>
#include <ctime>
#include <conio.h>

using namespace std;

void generator_test(double* point, double* center, int N, int n, double* center_test, int* n_point);
void identyfikacia(double* point, double* center, int N, int n, int* point_nr_center, double* max_dst);
void modyfikacja(double* point, double* center, int N, int n, int* point_nr_center);
void copy(int* point_nr_center_0, int* point_nr_center_1, int N);
int number_points_new_klaster(int* point_nr_center_0, int* point_nr_center_1, int N);
void means_distance(double* point, int N, int n, int* point_nr_center, double* means_dst);

double miara(double* x, double* c);

int main()
{
    const int N = 2000;  // number of points
    const int n = 5;     // number of centers 
    const int max_iter = 1000;
    double point[2][N];
    int point_nr_center_0[N];
    int point_nr_center_1[N];
    double center[2][n];

    double center_test[3][n] = { {0.6,0.7,0.9,0.2,0.2},
                                  {0.3,0.7,0.1,0.8,0.2},
                                  {0.1,0.3,0.1,0.2,0.15} };
    int nm_points[n] = { 400,500,300,300,500 };
    generator_test(*point, *center, N, n, *center_test, nm_points);

    ofstream file_0_points;
    file_0_points.open("points_test_v1.txt");
    for (int i = 0; i < N; i++)
        file_0_points << setw(12) << point[0][i] << "  " << setw(12) << point[1][i] << endl;
    file_0_points.close();

    ofstream file_0_centers;
    file_0_centers.open("centers_test_theor_v1.txt");
    for (int i = 0; i < n; i++)
    {
        file_0_centers << setw(10) << center_test[0][i] << "  " << setw(10) << center_test[1][i] << "  ";
        file_0_centers << setw(10) << center_test[2][i] << "  " << setw(10) << nm_points[i] << endl;
    }
    file_0_centers.close();

    ofstream file_1_centers;
    file_1_centers.open("centers_iteracie_v1.txt");

    double means_dst[n];
    for (int j = 0; j < n; j++)
        means_dst[j] = 1;
    identyfikacia(*point, *center, N, n, point_nr_center_0, means_dst);     // pochatkova inicializacja
    //    means_distance(*point,N,n,point_nr_center_0,means_dst);

    /////////////////////////////////////////////////////////////////////////////////////////////////////
    /////             iteracijnyj proces
    /////////////////////////////////////////////////////////////////////////////////////////////////////

    for (int iter = 0; iter < max_iter; iter++)
    {
        for (int j = 0; j < n; j++)
            file_1_centers << setw(12) << center[0][j] << "  " << setw(12) << center[1][j] << endl;
        file_1_centers << endl;
        modyfikacja(*point, *center, N, n, point_nr_center_0);
        identyfikacia(*point, *center, N, n, point_nr_center_1, means_dst);
        means_distance(*point, N, n, point_nr_center_1, means_dst);
        int k = number_points_new_klaster(point_nr_center_0, point_nr_center_1, N);
        cout << "nr_iter = " << setw(4) << iter << "       k= " << setw(6) << k << endl;
        if (k == 0) break;
        copy(point_nr_center_0, point_nr_center_1, N);
    }
    file_1_centers.close();

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////      vypysuvannia rezultativ
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    ofstream file_1_points_klasters;
    file_1_points_klasters.open("points_klasters_v1.txt");
    for (int k = 0; k < n; k++)
        for (int i = 0; i < N; i++)
            if (point_nr_center_1[i] == k)
            {
                file_1_points_klasters << setw(12) << point[0][i];
                file_1_points_klasters << setw(12) << point[1][i];
                file_1_points_klasters << setw(12) << point_nr_center_1[i] + 1 << endl;
            }
    file_1_points_klasters.close();

    return 0;
}

void identyfikacia(double* point, double* center, int N, int n, int* point_nr_center, double* means_dst)
{
    double* point_0 = new double[2];
    double* center_0 = new double[2];
    for (int i = 0; i < N; i++)
    {
        point_0[0] = point[0 * N + i];
        point_0[1] = point[1 * N + i];
        //center_0[0] = center[0];
        //center_0[1] = center[N];
        double miara_min = miara(point_0, center_0) / means_dst[0];
        for (int j = 0; j < n; j++)
        {
            center_0[0] = center[0 * n + j];
            center_0[1] = center[1 * n + j];
            //double miara_min;
            if (j == 0)
            {
                point_nr_center[i] = 0;
                miara_min = miara(point_0, center_0) / means_dst[j];
            }
            else
            {
                double temp = miara(point_0, center_0) / means_dst[j];
                if (miara_min > temp)
                {
                    miara_min = temp;
                    point_nr_center[i] = j;
                }
            }
        }
    }
}

void modyfikacja(double* point, double* center, int N, int n, int* point_nr_center)
{
    int* nr_point_center = new int[n];
    for (int j = 0; j < n; j++)
    {
        center[0 * n + j] = 0;
        center[1 * n + j] = 0;
        nr_point_center[j] = 0;
    }
    for (int i = 0; i < N; i++)
    {
        int temp = point_nr_center[i];
        center[0 * n + temp] += point[0 * N + i];
        center[1 * n + temp] += point[1 * N + i];
        nr_point_center[temp]++;
    }
    for (int j = 0; j < n; j++)
    {
        center[0 * n + j] = center[0 * n + j] / nr_point_center[j];
        center[1 * n + j] = center[1 * n + j] / nr_point_center[j];
    }
}

double miara(double* x, double* c)
{
    return sqrt(pow(x[0] - c[0], 2) + pow(x[1] - c[1], 2));
}

void copy(int* point_nr_center_0, int* point_nr_center_1, int N)
{
    for (int i = 0; i < N; i++)
        point_nr_center_0[i] = point_nr_center_1[i];
}

int number_points_new_klaster(int* point_nr_center_0, int* point_nr_center_1, int N)
{
    int k = 0;
    for (int i = 0; i < N; i++)
        if (point_nr_center_0[i] != point_nr_center_1[i]) k++;
    return k;
}

void generator_test(double* point, double* center, int N, int n, double* center_test, int* n_point)
{
    int* nr_k = new int[n];
    for (int k = 0; k < n; k++)
        nr_k[k] = 0;

    srand(131);
    int nr = 0;
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

    /*    for (int i=0;i<n;i++)
        for (int k=0;k<2;k++)
            center[k*n+i] = double(rand())/double(RAND_MAX);*/

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

void means_distance(double* point, int N, int n, int* point_nr_center, double* means_dst)
{
    int* nm = new int[n];
    for (int j = 0; j < n; j++)
    {
        means_dst[j] = 0;
        nm[j] = 0;
    }
    for (int j = 0; j < n; j++)
        for (int i = 0; i < N; i++)
            for (int k = i + 1; k < N; k++)
                if ((point_nr_center[i] == j) && (point_nr_center[k] == j))
                {
                    double P_1[2];
                    double P_2[2];
                    P_1[0] = point[0 * N + i];
                    P_1[1] = point[1 * N + i];
                    P_2[0] = point[0 * N + k];
                    P_2[1] = point[1 * N + k];
                    double temp = miara(P_1, P_2);
                    means_dst[j] += temp;
                    nm[j]++;
                }
    for (int j = 0; j < n; j++)
        means_dst[j] = means_dst[j] / nm[j];
}