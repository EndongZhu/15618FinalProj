#include <stdlib.h>
#include <time.h>
#include <float.h>
#include <mpi.h>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
using namespace std;

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


class Population {
public:
    int length;
    int capacity;
    int** individuals;
    int** next_generation;
    double mutate_rate;
    double** distance;
    int generation_num;

    Population(int capacity, int length, double** distance) {
        this->distance = distance;
        this->capacity = capacity;
        this->length = length;
        this->mutate_rate = 0.05;
        this->generation_num = 1;
        individuals = new int*[capacity];
        next_generation = new int*[capacity];
        for (int i = 0; i < capacity; i++) {
            next_generation[i] = new int[length];
            individuals[i] = new int[length];
            for (int j = 0; j < length; j++) {
                individuals[i][j] = j;
            }
            random_shuffle(&individuals[i][0], &individuals[i][length]);
            // for (int j = 0; j < length; j++) {
            //     cout << individuals[i][j] << " ";
            // }
            // cout << getIndividualCost(individuals[i]) << endl;
        }
        // cout << endl;
    };

    void crossover(int* parentA, int* parentB, int* child) {
        int start = rand() % length;
        int end = rand() % length;
        bool used[length];
        memset(used, 0, length*sizeof(bool));
        for (int i = start; i != end; i = (i+1)%capacity) {
            child[i] = parentA[i];
            used[parentA[i]] = 1;
        }
        int idx = end;
        for (int i = 0; i < length; i++) {
            if (used[parentB[i]] == 0) {
                child[idx] = parentB[i];
                used[parentB[i]] = 1;
                idx = (idx+1) % capacity;
            }
        }
    }

    void mutate(int* child) {
        for (int i = 0; i < length; i++) {
            double r = ((double) rand() / (RAND_MAX));
            if (r < mutate_rate) {
                int j = rand() % length;
                swap(child[i], child[j]);
            }
        }
    }

    void evolve() {
        // select the best individual as elite
        double best_score = DBL_MAX;
        for (int i = 0; i < capacity; i++) {
            double cur_cost = getIndividualCost(individuals[i]);
            if (cur_cost < best_score) {
                best_score = cur_cost;
                for (int j = 0; j < length; j++) {
                    next_generation[0][j] = individuals[i][j];
                }
            }
        }
        // using tournament select to generate children
        for (int i = 1; i < capacity; i++) {
            int a = tournamentSelect(capacity/8);
            int b = tournamentSelect(capacity/8);
            while (b == a) {
                b = tournamentSelect(capacity/8);
            }
            crossover(individuals[a], individuals[b], next_generation[i]);
            mutate(next_generation[i]);
        }
        // copy next_generation to population
        for (int i = 0; i < capacity; i++) {
            for (int j = 0; j < length; j++) {
                individuals[i][j] = next_generation[i][j];
                // if (i == 0)
                //     cout << next_generation[i][j] << " ";
            }
            // if (i == 0)
            //     cout << getIndividualCost(next_generation[i]) << endl;
        }
        // cout << endl;
        // mutate_rate decay
        generation_num += 1;
        mutate_rate *= (1 - 0.2 / generation_num);
        // cout << mutate_rate << endl;
    }

    int tournamentSelect(int tournamentNum) {
        int res = -1;
        double min_dist = DBL_MAX;
        for (int i = 0; i < tournamentNum; i++) {
            int idx = rand() % capacity;
            double cost = getIndividualCost(individuals[idx]);
            if (cost < min_dist) {
                res = idx;
                min_dist = cost;
            }
        }
        return res;
    }

    double getIndividualCost(int* tour) {
        double cost = 0;
        for (int i = 1; i < length; i++) {
            int a = tour[i];
            int b = tour[i-1];
            cost += distance[a][b];
        }
        cost += distance[tour[0]][tour[length-1]];
        return cost;
    }
};


class TSP {
public:
    int cityNum;
    City* cities;
    double** distance;
    Population* p;

    TSP(int num) {
        cityNum = num;
        cities = new City[num];
        distance = new double*[num];
        for (int i = 0; i < num; i++) {
            distance[i] = new double[num];
            for (int j = 0; j < num; j++) {
                distance[i][j] = calc_distance(cities[i], cities[j]);
                // cout << setprecision(4) << setw(5) << distance[i][j] << " ";
            }
            // cout << endl;
        }
        // cout << endl;
    };

    TSP(string filename) {
        ifstream infile(filename.c_str());
        infile >> cityNum;
        cities = new City[cityNum];
        int id;
        double x, y;
        while (infile >> id >> x >> y) {
            cities[id-1].x = x;
            cities[id-1].y = y;
        }
        distance = new double*[cityNum];
        for (int i = 0; i < cityNum; i++) {
            distance[i] = new double[cityNum];
            for (int j = 0; j < cityNum; j++) {
                distance[i][j] = calc_distance(cities[i], cities[j]);
            }
        }
        p = new Population(cityNum, cityNum, distance);
    }

private:
    double calc_distance(City a, City b) {
        double dist_x = (a.x - b.x) * (a.x - b.x);
        double dist_y = (a.y - b.y) * (a.y - b.y);
        return sqrt(dist_x + dist_y);
    }
};

int main(int argc, char *argv[]) {
    int process_count = 1;
    int this_zone = 0;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &process_count);
    MPI_Comm_rank(MPI_COMM_WORLD, &this_zone);

    string filename(argv[1]);
    srand(time(NULL) * this_zone);
    TSP t(filename);

    for (int round = 0; round < 10; round++) {
        for (int i = 0; i < 10000; i++) {
            t.p->evolve();
        }
        if (this_zone == 0) {
            int res[process_count][t.p->length];
            memcpy(res[0], t.p->individuals[0], sizeof(int) * t.p->length);
            for (int i = 1; i < process_count; i++) {
                MPI_Recv(res[i], t.p->length, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                memcpy(t.p->individuals[i], res[i], sizeof(int) * t.p->length);
            }
            MPI_Bcast(res, process_count * t.p->length, MPI_INT, 0, MPI_COMM_WORLD);
        } else {
            MPI_Send(t.p->individuals[0], t.p->length, MPI_INT, 0, 0, MPI_COMM_WORLD);
            int optimal[process_count * t.p->length];
            MPI_Bcast(optimal, process_count * t.p->length, MPI_INT, 0, MPI_COMM_WORLD);
            for (int i = 0; i < process_count; i++) {
                int start = i * t.p->length;
                for (int j = 0; j < t.p->length; j++) {
                    t.p->individuals[i][j] = optimal[start+j];
                }
            }
        }
    }

    for (int i = 0; i < 1000; i++) {
        t.p->evolve();
    }

    if (this_zone == 0) {
        cout << this_zone << " OPTIMAL:" ;
        for (int i = 0; i < t.p->length; i++)
            cout << t.p->individuals[0][i] << " " ;
        cout << t.p->getIndividualCost(t.p->individuals[0]) << endl;
    }

    MPI_Finalize();
    return 0;
}
