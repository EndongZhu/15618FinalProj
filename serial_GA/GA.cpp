#include <stdlib.h>
#include <time.h>
#include <cmath>
#include <algorithm>
#include <iostream>
using namespace std;

struct City {
    int x;
    int y;
    City() {
        x = rand() % 200;
        y = rand() % 200;
    }
    City(int new_x, int new_y) {
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

    Population(int capacity, int length, double** distance) {
        this->distance = distance;
        this->capacity = capacity;
        this->length = length;
        this->mutate_rate = 0.2;
        individuals = new int*[capacity];
        next_generation = new int*[capacity];
        for (int i = 0; i < capacity; i++) {
            next_generation[i] = new int[length];
            individuals[i] = new int[length];
            for (int j = 0; j < length; j++) {
                individuals[i][j] = j;
            }
            random_shuffle(&individuals[i][0], &individuals[i][length]);
            for (int j = 0; j < length; j++) {
                cout << individuals[i][j] << " ";
            }
            cout << getIndividualCost(individuals[i]) << endl;
        }
        cout << endl;
    };

    void crossover(int* parentA, int* parentB, int* child) {
        int start = rand() % length;
        int end = rand() % length;
        if (end < start) {
            swap(start, end);
        }
        bool used[length] = {0};
        for (int i = start; i <= end; i++) {
            child[i] = parentA[i];
            used[parentA[i]] = 1;
        }
        int idx = 0;
        for (int i = 0; i < length; i++) {
            if (idx == start) {
                idx = end+1;
            }
            if (used[parentB[i]] == 0) {
                child[idx++] = parentB[i];
                used[parentB[i]] == 1;
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
        for (int i = 0; i < capacity; i++) {
            int a = tournamentSelect(4);
            int b = tournamentSelect(4);
            crossover(individuals[a], individuals[b], next_generation[i]);
            mutate(next_generation[i]);
            for (int j = 0; j < length; j++) {
                cout << next_generation[i][j] << " ";
            }
            cout << getIndividualCost(next_generation[i]) << endl;
        }
        cout << endl;
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

private:
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

    TSP(int num) {
        cityNum = num;
        cities = new City[num];
        distance = new double*[num];
        for (int i = 0; i < num; i++) {
            distance[i] = new double[num];
            for (int j = 0; j < num; j++) {
                distance[i][j] = calc_distance(cities[i], cities[j]);
                cout << setprecision(4) << setw(5) << distance[i][j] << " ";
            }
            cout << endl;
        }
        cout << endl;
    };

private:
    double calc_distance(City a, City b) {
        double dist_x = (a.x - b.x) * (a.x - b.x);
        double dist_y = (a.y - b.y) * (a.y - b.y);
        return sqrt(dist_x + dist_y);
    }
};

int main() {
    TSP t(15);
    Population p(10, 15, t.distance);
    for (int i = 0; i < 5; i++) {
        p.evolve();
    }
}
