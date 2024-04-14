#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <cstdlib>
#include <ctime>
#include <conio.h>

using namespace std;

void generator_test(double** cities, int n);
double cities_distance(double** cities, int i, int j);
int ant_decision(double** cities, double* feromon_0, int* ant_kopia, int cities_0, int flag);
int random_decision(double* p, int n);
void min_2(int& i_0, int& j_0);
double min_vector(double* vector, int n);
double max_vector(double* vector, int n);
double means_vector(double* vector, int n);
void permutation(int* v, int n);
double f_st(double cykl);
int min_vector_int(int* vector, int n);

int nm_cities;
int nm_ant;
double alfa = 1;
double beta = 2;

int cykl_0;

int main()
{
    char  array_cities[40] = "array_cities.txt";
    char file_name_way[40] = "ants_way_cykl_30_s4.txt";
    char solution[40] = "SOLUTION_30_s4.txt";
    char graf_fig[40] = "graf_fig_30_s4.txt";
    int znak;

    double tau_0 = 1;
    double ro = 0.5;
    double tau = 0.1;

    int max_krok = 20000;
    int delta_s = 10;

    int nm_symulation;
    cout << "enter number symulation : ";
    cin >> nm_symulation;
    cout << endl << "number symulation is : " << setw(2) << nm_symulation << endl << endl;

    srand(9001);
    if (nm_symulation - 1)
    {
        cout << "enter number cities: ";
        cin >> nm_cities;
    }
    else
        nm_cities = 25 + int(11 * double(rand()) / double(RAND_MAX));
    cout << endl << "number cities = " << setw(2) << nm_cities << endl;

    double** cities;
    cities = new double* [nm_cities];
    for (int i_cities = 0; i_cities < nm_cities; i_cities++)
        cities[i_cities] = new double[2];           // masyw koordynat mist

    if (nm_symulation - 1)
    {
        ifstream file_read;
        file_read.open(array_cities);                   // koordynaty mist do fajlu
        if (!file_read.good())
        {
            cout << "ERROR READ: FILE READ NO";
            cout << "ENTER KEYBOARD KEA 1 AND ENTER" << endl;
            cin >> znak;
            return 1;
        }
        for (int i_cities = 0; i_cities < nm_cities; i_cities++)
            file_read >> cities[i_cities][0] >> cities[i_cities][1];
        file_read.close();
    }
    else
    {

        generator_test(cities, nm_cities);               // stvoriujemo mista
        ofstream file_cities;
        file_cities.open(array_cities);                   // koordynaty mist do fajlu
        for (int i_cities = 0; i_cities < nm_cities; i_cities++)
            file_cities << setw(10) << cities[i_cities][0] << "   " << setw(10) << cities[i_cities][1] << endl;
        file_cities.close();
    }

    if (nm_symulation - 1)
    {
        for (int i_cities = 0; i_cities < nm_cities; i_cities++)
        {
            for (int j_cities = i_cities + 1; j_cities < nm_cities; j_cities++)
            {
                int temp = int(cities_distance(cities, i_cities, j_cities));
                if ((temp < 10) || (temp > 100))
                {
                    cout << "ERROR READ: DISTANCE : " << setw(4) << temp << endl;
                    cin >> znak;
                    cout << "ENTER KEYBOARD KEA 1 AND ENTER" << endl;
                    return 1;
                }
                else temp = 0;
            }
        }
    }

    cout << "enter number ant (<= number cities): ";
    cin >> nm_ant;
    cout << endl << "number ant = " << setw(2) << nm_ant << endl << endl;

    double* feromon_0;
    feromon_0 = new double[nm_cities * (nm_cities - 1) / 2];          // masyv feromonu na stezhci
    double* feromon_1;
    feromon_1 = new double[nm_cities * (nm_cities - 1) / 2];           // masyv feromonu v ostanniomu cykli

    int* nm_ant_way;
    nm_ant_way = new int[nm_cities * (nm_cities - 1) / 2];             // masyv kilkosti murach, szczo vidvidaly stezhku
    int* ant_flag;
    ant_flag = new int[nm_ant];

    for (int i_pozycia = 0; i_pozycia < nm_cities * (nm_cities - 1) / 2; i_pozycia++)
    {
        feromon_0[i_pozycia] = tau_0;                            // pochatkovyj feromon                           
        feromon_1[i_pozycia] = 0;
        nm_ant_way[i_pozycia] = 0;
    }
    for (int i_ant = 0; i_ant < nm_ant; i_ant++)
        ant_flag[i_ant] = 0;

    ////////////////////////////////////////////////////////////////////////////////////////////////
    //////            WYPYSUVANNIA PRY TESTUVANNI ROBOTY PROHRAMY
    ////////////////////////////////////////////////////////////////////////////////////////////////    
    if (nm_cities < 11)                                       // testuvannia: vidstani mizh mistamy na ekran
    {
        cout << "CITIES DISTANCE" << endl << endl;
        for (int i_cities = 0; i_cities < nm_cities; i_cities++)
        {
            for (int j_cities = 0; j_cities < nm_cities; j_cities++)
            {
                int temp;
                if (i_cities != j_cities) temp = int(cities_distance(cities, i_cities, j_cities));
                else temp = 0;
                cout << setw(3) << temp << "  ";
            }
            cout << endl;
        }
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    int** ant_way_0;                          // masyv z informacijeju z ostannioho cyklu
    int** ant_way_1;                          // bizhuczyj masyv mist
    int** ant_way_2;                          // masyv blokuvannia + bizhucha informacia v cykli
    ant_way_0 = new int* [nm_ant];
    ant_way_1 = new int* [nm_ant];
    ant_way_2 = new int* [nm_ant];

    for (int i_ant = 0; i_ant < nm_ant; i_ant++)
    {
        ant_way_0[i_ant] = new int[nm_cities + 3];
        ant_way_1[i_ant] = new int[nm_cities];
        ant_way_2[i_ant] = new int[nm_cities + 3];
    }

    int* ant_cities_poch = new int[nm_ant];
    if (nm_ant < nm_cities)
    {
        int* temp_poch = new int[nm_cities];
        for (int i_c = 0; i_c < nm_cities; i_c++)
            temp_poch[i_c] = i_c;
        permutation(temp_poch, nm_cities);
        for (int i_ant = 0; i_ant < nm_ant; i_ant++)
            ant_cities_poch[i_ant] = temp_poch[i_ant];
    }
    else if (nm_ant == nm_cities)
    {
        for (int i_ant = 0; i_ant < nm_ant; i_ant++)
            ant_cities_poch[i_ant] = i_ant;
    }
    else
    {
        cout << "ERROR: NM_ANT > NM_CITIES" << endl;
        //getch();
        return 1;
    }

    for (int i_ant = 0; i_ant < nm_ant; i_ant++)
    {
        ant_way_0[i_ant][nm_cities + 0] = 0;          //  projdenyj v cykli szlach 
        ant_way_0[i_ant][nm_cities + 1] = 0;          //  ves szlach
        ant_way_0[i_ant][nm_cities + 2] = 0;          //  kilkist' projdenych cykliv
        ant_way_2[i_ant][nm_cities + 0] = 0;          //  avans odynyc szlachu
        ant_way_2[i_ant][nm_cities + 1] = 0;          //  skilky projdeno mist v cykli (narazi 0)
        ant_way_2[i_ant][nm_cities + 2] = 0;          //  bizhuchyj projdenyj szlach v cykli

        for (int i_cities = 0; i_cities < nm_cities; i_cities++)
        {
            if (i_cities == 0)      ant_way_1[i_ant][i_cities] = ant_cities_poch[i_ant];      // misto startu murachy
            if (i_cities != ant_cities_poch[i_ant]) ant_way_2[i_ant][i_cities] = 1;          // inszi mista vilni
            else                  ant_way_2[i_ant][i_cities] = 0;          // misto startu blokujemo
        }
    }
    int* ant_kopia = new int[nm_cities + 3];          // masyv [] dla kopii riadka masyvu ant_way_2[][]
    int* ant_cykl = new int[nm_ant];
    
    for (int i_ant = 0; i_ant < nm_ant; i_ant++)             // stan murachy v cykli 0 (tryvaje) 1 (muracha projszla cykl) 
        ant_cykl[i_ant] = 0;

    cykl_0 = 0;                                      // nomer cyklu
    double cykl_flag = 0;                            // stan cyklu (pochynajet'sia - tryvaje)
    double minw_way = 0;                             // najkorotszyj szlach v zahalnomu cykli
    double means_way = 0;                            // serednij szlach v zahalnomu cykli
    double maxw_way = 0;                             // najkorotszyj szlach v zahalnomu cykli
    int temp;
    int nm_ant_way_no_zero;
    ofstream file_way;
    file_way.open(file_name_way);
    int cykl_st = 160;

    for (int i_krok = 1; i_krok <= max_krok; i_krok++)      // start zahalnyj
    {
        for (int i_ant = 0; i_ant < nm_ant; i_ant++)        // start murachy z nomerom i_ant 
        {
            ant_way_2[i_ant][nm_cities + 0] += delta_s;  // avans odynyc szlachu >0 vyruszjemo z mista
            // <0 narazi v dorozi do nastupnoho mista
            if (ant_way_2[i_ant][nm_cities + 0] > 0)       // mozhemo ruszaty z nastupnoho mista
            {
                int ant_cities_nr = ant_way_2[i_ant][nm_cities + 1];     // fiksujemo z jakoho po nomeru mista vychodymo
                int ant_cities_0 = ant_way_1[i_ant][ant_cities_nr];   // nomer mista z jakoho vychodymo  
                int ant_cities_1;                                      // zminna - nomer mista do jakoho pidemo

                if (ant_cities_nr < nm_cities - 2)             // lyszylos' projty 2 abo bilsze mist
                {
                    for (int i = 0; i < nm_cities + 3; i++)
                        ant_kopia[i] = ant_way_2[i_ant][i];     // robymo kopiju masywu blokuvannia - dla peredach do 
                    // do funkcii pryniattia riszennia
                    ant_cities_1 = ant_decision(cities, feromon_0, ant_kopia, ant_cities_0, ant_flag[i_ant]);  // pryjmajemo riszennia
                    ant_way_2[i_ant][nm_cities + 1]++;            // zbilszujemo chyslo projdenych mist
                    ant_cities_nr = ant_way_2[i_ant][nm_cities + 1];
                    ant_way_1[i_ant][ant_cities_nr] = ant_cities_1;  // dodajemo misto do bizhuchoho spysku
                    ant_way_2[i_ant][ant_cities_1] = 0;             // blokujemo na majbutnie
                }
                else if (ant_cities_nr == nm_cities - 2)             // zalyszylos' tilky odne nevidvidane misto
                {
                    for (int i_cities = 0; i_cities < nm_cities; i_cities++)     // szukajemo ostannie nezablokovane
                        if (ant_way_2[i_ant][i_cities])
                            ant_cities_1 = i_cities;                   // jdemo do ostannioho nezablokovanoho
                    ant_way_2[i_ant][nm_cities + 1]++;                       // zbilszujemo chyslo prpjdenych mist
                    ant_way_1[i_ant][nm_cities - 1] = ant_cities_1;          // dodajemo misto do bizhuchoho spysku
                }
                else                        // vsi mista projdeno - vertajemos' do pochatkovoho mista cyklu
                {
                    temp = int(cities_distance(cities, ant_cities_0, ant_cities_poch[i_ant]));   // 
                    // temp - vidstan', jaku treba projty
                    ant_cities_1 = ant_cities_poch[i_ant];                                     // jdemo do pochatkovoho mista
                    for (int i_cities = 0; i_cities < nm_cities; i_cities++)
                    {                                 // povernulys' do startu: kopijujemo poslidovnist' projdenych mist
                                                      // nanovo rozblokovujemo mista dla nastupnoho cyklu
                        ant_way_0[i_ant][i_cities] = ant_way_1[i_ant][i_cities];
                        if (ant_cities_poch[i_ant] != i_cities)  ant_way_2[i_ant][i_cities] = 1;
                        else                  ant_way_2[i_ant][i_cities] = 0;
                    }
                    // slach v cykli : do bizhuchoho szlachu + ostannij promizhok slachu v cykli 
                    ant_way_0[i_ant][nm_cities + 0] = ant_way_2[i_ant][nm_cities + 2] + temp;   // projdenyj v cykli szlach 
                    ant_way_0[i_ant][nm_cities + 1] += ant_way_0[i_ant][nm_cities + 0];  // slach u vsich cyklach
                    ant_way_0[i_ant][nm_cities + 2]++;                                 // kilkist' cykliv
                    ant_way_2[i_ant][nm_cities + 2] = -temp;   // zanulujemo bizhuchyj projdenyj szlach
                    // pisla else{...} dodamo temp : tomu vyjde 0
                    ant_way_2[i_ant][nm_cities + 1] = 0;       // v nastupnomu cykli projdeno 0 mist
                    // murachy svij cykl zakinchujut' v riznyj moment
                    // povnyj cykl zakinchyt'sia, koly ostannia muracha zakinchyt' cykl
                    ant_cykl[i_ant] = 1;   // muracha z nomerom i_ant zakinchyla cykl
                }
                temp = int(cities_distance(cities, ant_cities_0, ant_cities_1));   // projdenyj szlach vid ant_cities_0 do ant_cities_1
                ant_way_2[i_ant][nm_cities + 0] -= temp;     // kilkist odynyc' szlachu   
                ant_way_2[i_ant][nm_cities + 2] += temp;     // szukajemo bizhuczyj projdenyj szlach v cykli
                min_2(ant_cities_0, ant_cities_1);         // szukajemo nomer stezhky
                temp = (2 * nm_cities - ant_cities_0 - 1) * ant_cities_0 / 2 + ant_cities_1 - ant_cities_0 - 1;
                nm_ant_way[temp]++;    // muracha vidvidala stezku : fiksujemo v masyvi nm_ant_way
                if (cykl_0 > 0)
                    feromon_1[temp] += (tau + 0.5 * (1 - exp(-0.001 * cykl_0))) * pow(means_way / double(ant_way_0[i_ant][nm_cities]), 3);
                else
                    feromon_1[temp] += tau;
            }
        }

        int nm_ant_cykl = 0;   // skilky murach zakinchyly bizhuchyj cykl
        for (int i_a = 0; i_a < nm_ant; i_a++)
            nm_ant_cykl += ant_cykl[i_a];
        if (nm_ant_cykl == nm_ant)   // vsi murachy majut' zakinchenyj cykl : perechodymo do novoho
        {
            cykl_flag = 1;
            cykl_0++;
            for (int i_a = 0; i_a < nm_ant; i_a++)
                ant_cykl[i_a] = 0;
        }

        if (cykl_flag)
        {
            cykl_flag = 0; // vychodymo z if {...} i pochynajemo novyj cykl
            nm_ant_way_no_zero = 0;
            for (int i = 0; i < nm_cities * (nm_cities - 1) / 2; i++)
                if (nm_ant_way[i] != 0) nm_ant_way_no_zero++;    // skilky riznych stezhok vidvidano v cykli

            double* ant_way_00 = new double[nm_ant];   // fiksujemo szlach v cykli kozhnoji murachy  
            for (int i = 0; i < nm_ant; i++)
                ant_way_00[i] = double(ant_way_0[i][nm_cities]);
            minw_way = min_vector(ant_way_00, nm_ant);    // najkorotszyj szlach : vykorystajemo z perszoho cyklu
            means_way = means_vector(ant_way_00, nm_ant); // serednij
            maxw_way = max_vector(ant_way_00, nm_ant);    // najdovszyj

            // dla obczyslennia feromonu : v majbutniomu dovszyj szlach
            // mensze feromonu

            for (int i_ant0 = 0; i_ant0 < nm_ant; i_ant0++)
            {
                int temp_1 = ant_way_0[i_ant0][nm_cities] < means_way; // + 0.5*(maxw_way-means_way);
                int temp_2 = (maxw_way - minw_way) / maxw_way < 0.1;
                if (temp_1 || temp_2)
                    ant_flag[i_ant0] = 0;
                else
                    ant_flag[i_ant0] = 1;
                if (ant_way_0[i_ant0][nm_cities] == maxw_way) ant_flag[i_ant0] = 1;
            }

            if (cykl_0 > int(cykl_st) + nm_ant)
                for (int i_ant0 = 0; i_ant0 < nm_ant; i_ant0++)
                    if (ant_way_0[i_ant0][nm_cities] == minw_way)
                    {
                        int cit_0;
                        int cit_1;
                        int i_poz;
                        for (int i_cit = 0; i_cit < nm_cities - 1; i_cit++)
                        {
                            cit_0 = ant_way_0[i_ant0][i_cit];
                            cit_1 = ant_way_0[i_ant0][i_cit + 1];
                            min_2(cit_0, cit_1);         // szukajemo nomer stezhky
                            i_poz = (2 * nm_cities - cit_0 - 1) * cit_0 / 2 + cit_1 - cit_0 - 1;
                            feromon_0[i_poz] = f_st(cykl_0);
                        }
                        cit_0 = ant_way_0[i_ant0][nm_cities - 1];
                        cit_1 = ant_way_0[i_ant0][0];
                        min_2(cit_0, cit_1);         // szukajemo nomer stezhky
                        i_poz = (2 * nm_cities - cit_0 - 1) * cit_0 / 2 + cit_1 - cit_0 - 1;
                        feromon_0[i_poz] = f_st(cykl_0);
                    }

            double* fer_temp = new double[nm_cities * (nm_cities - 1) / 2];
            for (int i = 0; i < nm_cities * (nm_cities - 1) / 2; i++)
                if (feromon_0[i] < 0.9 * f_st(cykl_st))
                    fer_temp[i] = feromon_0[i];
                else
                    fer_temp[i] = 0.1;
            double feromon_m = ro * means_vector(fer_temp, nm_cities * (nm_cities - 1) / 2);
            for (int i = 0; i < nm_cities * (nm_cities - 1) / 2; i++)   // obnovlujemo feromon
            {
                feromon_0[i] -= feromon_m;
                if (feromon_0[i] < 0) feromon_0[i] = 0;
                feromon_0[i] += feromon_1[i];
                feromon_1[i] = 0;
            }
            for (int i_pozycia = 0; i_pozycia < nm_cities * (nm_cities - 1) / 2; i_pozycia++)
                nm_ant_way[i_pozycia] = 0;

            if (!((cykl_0 - 1) % 25))
            {
                file_way << "cykl = " << setw(5) << cykl_0 << endl;
                cout << "cykl = " << setw(5) << cykl_0 << endl;

                for (int i_ant = 0; i_ant < nm_ant; i_ant++)
                {
                    file_way << " " << setw(4) << ant_way_0[i_ant][nm_cities];  // do fajlu i na ekran szlach murachy 
                    cout << " " << setw(4) << ant_way_0[i_ant][nm_cities];    // v ostatnio projdenomu cykli                                                             
                    if (!((i_ant + 1) % 10)) file_way << endl;
                    if (!((i_ant + 1) % 10)) cout << endl;
                }
                file_way << endl << endl;
                cout << endl << endl;
                //                   getch();
            }
            if (nm_ant_way_no_zero == nm_cities) break;
        }
    }

    file_way << "cykl = " << setw(5) << cykl_0 << endl;
    cout << "cykl = " << setw(5) << cykl_0 << endl;
    for (int i_ant = 0; i_ant < nm_ant; i_ant++)
    {
        file_way << " " << setw(4) << ant_way_0[i_ant][nm_cities];  // do fajlu szlach murachy 
        cout << " " << setw(4) << ant_way_0[i_ant][nm_cities];  // v ostatnimu cykli
        if (!((i_ant + 1) % 10)) cout << endl;
        if (!((i_ant + 1) % 10)) file_way << endl;
    }
    file_way << endl << endl;
    cout << endl << endl;
    file_way.close();

    ofstream file_solution;
    ofstream file_fig;
    //    file_solution.open(file_solution,ios::app|ios::out);
    file_solution.open(solution);
    file_fig.open(graf_fig);
    for (int i_ant = 0; i_ant < nm_ant; i_ant++)
        ant_kopia[i_ant] = ant_way_0[i_ant][nm_cities];
    int min_way = min_vector_int(ant_kopia, nm_ant);

    cykl_flag = 1;
    for (int i_ant = 0; i_ant < nm_ant; i_ant++)
    {
        if ((ant_way_0[i_ant][nm_cities] == min_way) && cykl_flag)
        {
            cout << "SOLUTION_WAY = " << setw(4) << ant_way_0[i_ant][nm_cities] << endl << endl;
            cout << "WAY: " << endl;
            file_solution << "NUMBER CITIES = " << setw(4) << nm_cities << endl;
            file_solution << "NUMBER ANT    = " << setw(4) << nm_ant << endl << endl;
            file_solution << "SOLUTION_WAY =  " << setw(4) << ant_way_0[i_ant][nm_cities] << endl << endl;

            file_solution << "WAY: " << endl;
            int cities_0;
            int cities_1;
            int way_c0_c1;
            int way_kontrol = 0;
            for (int i_cities = 0; i_cities < nm_cities - 1; i_cities++)
            {
                cities_0 = ant_way_0[i_ant][i_cities];
                cities_1 = ant_way_0[i_ant][i_cities + 1];
                way_c0_c1 = int(cities_distance(cities, cities_0, cities_1));
                way_kontrol += way_c0_c1;
                cout << "WAY NR: " << setw(2) << i_cities + 1;
                cout << "  CITIES_0 NR: " << setw(2) << ant_way_0[i_ant][i_cities];
                cout << "  CITIES_1 NR: " << setw(2) << ant_way_0[i_ant][i_cities + 1];
                cout << "  L: " << setw(4) << way_c0_c1 << endl;
                file_solution << "WAY NR: " << setw(2) << i_cities + 1;
                file_solution << "  CITIES_0 NR: " << setw(2) << ant_way_0[i_ant][i_cities];
                file_solution << "  CITIES_1 NR: " << setw(2) << ant_way_0[i_ant][i_cities + 1];
                file_solution << "  L: " << setw(4) << way_c0_c1 << endl;
                file_fig << setw(10) << cities[cities_0][0] << "  ";
                file_fig << setw(10) << cities[cities_0][1] << endl;
            }
            cities_0 = ant_way_0[i_ant][nm_cities - 1];
            cities_1 = ant_way_0[i_ant][0];
            file_fig << setw(10) << cities[cities_0][0] << "  ";
            file_fig << setw(10) << cities[cities_0][1] << endl;
            file_fig << setw(10) << cities[cities_1][0] << "  ";
            file_fig << setw(10) << cities[cities_1][1] << endl;
            way_c0_c1 = int(cities_distance(cities, cities_0, cities_1));
            way_kontrol += way_c0_c1;
            cout << "WAY NR: " << setw(2) << nm_cities;
            cout << "  CITIES_0 NR: " << setw(2) << ant_way_0[i_ant][nm_cities - 1];
            cout << "  CITIES_1 NR: " << setw(2) << ant_way_0[i_ant][0];
            cout << "  L: " << setw(4) << way_c0_c1 << endl;
            cout << "                                      SUMA_WAY : " << setw(4) << way_kontrol << endl;
            file_solution << "WAY NR: " << setw(2) << nm_cities;
            file_solution << "  CITIES_0 NR: " << setw(2) << ant_way_0[i_ant][nm_cities - 1];
            file_solution << "  CITIES_1 NR: " << setw(2) << ant_way_0[i_ant][0];
            file_solution << "  L: " << setw(4) << way_c0_c1 << endl;
            file_solution << "                                      SUMA_WAY : " << setw(4) << way_kontrol << endl;
            cykl_flag = 0;
        }
    }
    file_solution << endl << endl;
    file_solution.close();

   // getch();
    return 0;
}

void min_2(int& i_0, int& j_0)
{
    if (i_0 > j_0)
    {
        int temp = i_0;
        i_0 = j_0;
        j_0 = temp;
    }
}

double cities_distance(double** cities, int i, int j)
{
    double temp = sqrt(pow(cities[i][0] - cities[j][0], 2) + pow(cities[i][1] - cities[j][1], 2));
    return temp;
}

int random_decision(double* p, int n)
{
    //    srand(time(NULL));
    srand(1001);
    double temp = double(rand()) / double(RAND_MAX);
    double temp_p = p[0] + 0.00000001;
    int decision = n - 1;
    for (int i = 0; i < n - 1; i++)
    {
        if (temp < temp_p)
        {
            decision = i;
            break;
        }
        temp_p += p[i + 1];
    }
    return decision;
}

int ant_decision(double** cities, double* feromon_0, int* ant_kopia, int cities_0, int flag)
{
    int pozycia = nm_cities - 1 - ant_kopia[nm_cities + 1];
    double* ant_fi = new double[pozycia];
    int* ant_cities = new int[pozycia];
    double ant_fi_suma = 0;
    int nr = 0;
    for (int i_cities = 0; i_cities < nm_cities; i_cities++)
        if (ant_kopia[i_cities])
        {
            int cities_01 = i_cities;
            int cities_00 = cities_0;
            min_2(cities_00, cities_01);
            int pozycia_0 = (2 * nm_cities - cities_00 - 1) * cities_00 / 2 + cities_01 - cities_00 - 1;
            double temp_1 = feromon_0[pozycia_0];
            double temp_2 = 10 / cities_distance(cities, cities_0, i_cities);
            ant_fi[nr] = exp(alfa * log(temp_1)) * exp(beta * log(temp_2));
            /*             if ((cykl_0==100)||(cykl_0==101)||(cykl_0==102))
                         {
                             cout << "i_c= " << setw(2) << i_cities << "  t_1= " << setw(10) << temp_1 << "  ";
                             cout << "  t_2= " << setw(10) << temp_2 << "  " << setw(10) << "fi= " << setw(10) << ant_fi[nr] << endl;
                         }*/
            ant_fi_suma += ant_fi[nr];
            ant_cities[nr] = i_cities;
            nr++;
        }
    for (int i_nr = 0; i_nr < nr; i_nr++)
        ant_fi[i_nr] = ant_fi[i_nr] / ant_fi_suma;
    int ant_cities_1;
    if (!flag)
    {
        double max_ant_fi = ant_fi[0];
        int max_inr = 0;
        for (int i_nr = 1; i_nr < nr; i_nr++)
            if (max_ant_fi < ant_fi[i_nr])
            {
                max_ant_fi = ant_fi[i_nr];
                max_inr = i_nr;
            }
        ant_cities_1 = ant_cities[max_inr];
    }
    else
        ant_cities_1 = ant_cities[random_decision(ant_fi, nr)];
    delete[] ant_fi;
    delete[] ant_cities;
    /*     if ((cykl_0==100)||(cykl_0==101)||(cykl_0==102))
         {
             cout << "decyzia : z  " << setw(2) << cities_0;
             cout << "  do   " << setw(2) << ant_cities_1 << endl << endl << endl;
             getch();
         }*/
    return ant_cities_1;
}

int miara(double* p1, double* p2)
{
    int temp = int(sqrt(pow(p1[0] - p2[0], 2) + pow(p1[1] - p2[1], 2)));
    return temp;
}

void generator_test(double** cities, int n)
{
    srand(131);
    int nr = 0;
    double p1[2];
    do
    {
        p1[0] = 100 * double(rand()) / double(RAND_MAX);
        p1[1] = 100 * double(rand()) / double(RAND_MAX);
        double temp = pow(p1[0] - 50, 2) + pow(p1[1] - 50, 2);
        if (temp < 2500)
        {
            int i_nr;
            if (nr == 0)
            {
                cities[nr][0] = p1[0];
                cities[nr][1] = p1[1];
                nr++;
            }
            else
            {
                for (i_nr = 0; i_nr < nr; i_nr++)
                {
                    double p2[2];
                    p2[0] = cities[i_nr][0];
                    p2[1] = cities[i_nr][1];
                    double temp = miara(p1, p2);
                    if ((temp < 10) || (temp > 100)) break;
                }
                if (i_nr == nr)
                {
                    cities[nr][0] = p1[0];
                    cities[nr][1] = p1[1];
                    nr++;
                }
            }
        }
    } while (nr < n);
}

double min_vector(double* vector, int n)
{
    double ekstr = vector[0];
    if (n == 0) return ekstr;
    for (int i = 1; i < n; i++)
        if (ekstr > vector[i])
            ekstr = vector[i];
    return ekstr;
}

int min_vector_int(int* vector, int n)
{
    int ekstr = vector[0];
    if (n == 0) return ekstr;
    for (int i = 1; i < n; i++)
        if (ekstr > vector[i])
            ekstr = vector[i];
    return ekstr;
}

double max_vector(double* vector, int n)
{
    double ekstr = vector[0];
    if (n == 0) return ekstr;
    for (int i = 1; i < n; i++)
        if (ekstr < vector[i])
            ekstr = vector[i];
    return ekstr;
}

double means_vector(double* vector, int n)
{
    double s = 0;
    for (int i = 0; i < n; i++)
        s += vector[i];
    return s / double(n);
}

void permutation(int* v, int n)
{
    int* v_0 = new int[n];
    for (int i = 0; i < n; i++)
        v_0[i] = 1;
    int i_nr = 0;
    do
    {
        int temp = rand() % n;
        if (v_0[temp])
        {
            v[i_nr] = temp;
            v_0[temp] = 0;
            i_nr++;
        }
    } while (i_nr < n);

    for (int i = 0; i < n; i++)
        if (v_0[i])
            v[n - 1] = i;
    delete[] v_0;
}

double f_st(double cykl)
{
    return exp(log(cykl / 16.) * 5);
}