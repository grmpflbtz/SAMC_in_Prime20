#ifndef RND_HPP
#define RND_HPP

#include <random>

// Declare engine - single instance for the whole code
extern std::mt19937 my_rng;

long int Seed(long int seed);
double RND();

#endif RND_HPP