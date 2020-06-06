#pragma once

#include <vector>
#include <random>

// generate k numbers from [a, b] - closed interval
// a should be non-negative, usually 0 or 1
std::vector<int> combination(const int &a, const int &b, const int &k, std::mt19937& rng);

struct uid_wrapper {
	//*
	int from, len;
	std::mt19937& rng;
	unsigned completePart;
	/*/
	std::mt19937& rng;
	std::uniform_int_distribution<int> uid;
	//*/
	uid_wrapper(int, int, std::mt19937&);

	int operator()();
};