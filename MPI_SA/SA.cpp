#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>
#include<math.h>
#include <mpi.h>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;

#define T0 10000.0  // Initial temperature
#define T_end 0.2
#define q  0.99   // Annealing factor
#define L 100  // Iterations at each temperature


struct City {
    double x;
    double y;
    City() {
        x = rand() % 200;
        y = rand() % 200;
    }
    City(double new_x, double new_y) {
        x = new_x;
        y = new_y;
    }
};


class SA {
public:
    double** distance;
    int city_num;
    City* cities;
    int* solution;
    int* new_solution;

    SA(int num) {
        city_num = num;
        cities = new City[num];
        distance = new double*[num];
        solution  = new int[num];
        new_solution  = new int[num];
        for (int i = 0; i < num; i++) {
            distance[i] = new double[num];
            for (int j = 0; j < num; j++) {
                distance[i][j] = calc_distance(cities[i], cities[j]);
            }
        }
        for (int i = 0; i < num; i++) {
            solution[i] = i;
        }
        random_shuffle(solution, solution + city_num);
    }


    SA(string filename) {
        ifstream infile(filename.c_str());
        infile >> city_num;
        cities = new City[city_num];
        solution  = new int[city_num];
        new_solution  = new int[city_num];
        int id;
        double x, y;
        while (infile >> id >> x >> y) {
            cities[id-1].x = x;
            cities[id-1].y = y;
        }
        distance = new double*[city_num];
        for (int i = 0; i < city_num; i++) {
            distance[i] = new double[city_num];
            for (int j = 0; j < city_num; j++) {
                distance[i][j] = calc_distance(cities[i], cities[j]);
            }
        }
        for (int i = 0; i < city_num; i++) {
            solution[i] = i;
        }
        random_shuffle(solution, solution + city_num);
    }


    double calc_distance(City a, City b) {
        double dist_x = (a.x - b.x) * (a.x - b.x);
        double dist_y = (a.y - b.y) * (a.y - b.y);
        return sqrt(dist_x + dist_y);
    }

    double get_path_len(int* tour) {
        double cost = 0;
        for (int i = 1; i < city_num; i++) {
            int a = tour[i];
            int b = tour[i-1];
            cost += distance[a][b];
        }
        cost += distance[tour[0]][tour[city_num-1]];
        return cost;
    }

    // Generate new solution
    void create_new() {
        memcpy(new_solution, solution, city_num * sizeof(int));
        double r1 = ((double)rand())/(RAND_MAX+1.0);
        double r2 = ((double)rand())/(RAND_MAX+1.0);
        // two positions to be swapped
        int pos1 = (int)(city_num*r1);
        int pos2 = (int)(city_num*r2);
        while(pos1 == pos2) {
            r2 = ((double)rand())/(RAND_MAX+1.0);
            pos2 = (int)(city_num*r2);
        }
        swap(new_solution[pos1], new_solution[pos2]);
    }

    void accept_new() {
        memcpy(solution, new_solution, city_num * sizeof(int));
    }

    void print_solution() {
        for (int i = 0; i < city_num; i++) {
            cout << solution[i] << " ";
        }
        cout << get_path_len(solution) << endl;
    }
};


int main(int argc, char **argv)
{
    int process_count = 1;
    int this_zone = 0;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &process_count);
    MPI_Comm_rank(MPI_COMM_WORLD, &this_zone);

    string filename(argv[1]);
    srand(time(NULL) * this_zone); //initialize random seed
    double T = T0;
    int count = 0; // record annealing numbers
    SA s(filename);

    while(T > T_end)
    {
        for(int i = 0; i < L; i++)
        {
            s.create_new(); // generate new solution
            double f1 = s.get_path_len(s.solution);
            double f2 = s.get_path_len(s.new_solution);
            double df = f2 - f1;
            // Metropolis fomula
            if(df > 0)
            {
                double r = ((double) rand() / (RAND_MAX));
                if(r <= exp(-df / T)) // accept new solution
                {
                    s.accept_new();
                }
            } else {
                s.accept_new();
            }
            count++;
        }
        T *= q; // annealing

        if (count >= 100) {
            if (this_zone == 0) {
                int res[process_count][s.city_num];
                memcpy(res[0], s.solution, sizeof(int) * s.city_num);
                double best_cost = s.get_path_len(res[0]);
                int best_idx = 0;
                for (int i = 1; i < process_count; i++) {
                    MPI_Recv(res[i], s.city_num, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    double cost = s.get_path_len(res[i]);
                    if (cost < best_cost) {
                        best_cost = cost;
                        best_idx = i;
                    }
                }
                MPI_Bcast(res[best_idx], s.city_num, MPI_INT, 0, MPI_COMM_WORLD);
            } else {
                MPI_Send(s.solution, s.city_num, MPI_INT, 0, 0, MPI_COMM_WORLD);
                MPI_Bcast(s.solution, s.city_num, MPI_INT, 0, MPI_COMM_WORLD);
            }
            count = 0;
        }
    }



    cout << "Process " << this_zone << ": ";
    s.print_solution();

    MPI_Finalize();
    return 0;
}
