#include <stdio.h>
#include <iostream>
#include <chrono>
#include "random.hpp"

std::mt19937 my_rng {};     // defines an engine

// Function to seed the random number generator from main file
// Negative value uses TIME as random seed
long unsigned int Seed(int seed) {
    if (seed < 0) {
        long rseed=static_cast<long unsigned int>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
        //std::cerr << "Randomizing random generator, seed is "<<rseed<<std::endl;
        my_rng.seed(rseed - seed);
        return rseed - seed;
    } else {
        //std::cerr << "User-provided seed is "<<seed<<std::endl;
        my_rng.seed(seed);
        return seed;
    }
}
double RND() {
    return my_rng();
}