#ifndef METAGEM_H
#define METAGEM_H

#include <iostream>
#include <fstream>
#include <chrono>
#include <unordered_map>
#include <unordered_set>

#include <boost/math/distributions/chi_squared.hpp>

#include "../include/CLI11-2.2.0/CLI11.hpp"
#include "../include/sparsepp-2018.2/spp.h"
#include "matrix_utils.h"
#include "read_parameters.h"
#include "file.h"

#define VERSION "1.0"

using std::cout;
using std::endl;
using std::cerr;
using spp::sparse_hash_map;


void metagem(CommandLine cmd);

template <class T>
inline void hash_combine(std::size_t & seed, const T & v)
{
  std::hash<T> hasher;
  seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

namespace std
{
  template<typename S, typename T> struct hash<pair<S, T>>
  {
    inline size_t operator()(const pair<S, T> & v) const
    {
      size_t seed = 0;
      ::hash_combine(seed, v.first);
      ::hash_combine(seed, v.second);
      return seed;
    }
  };
}

void printWelcome();
void printMetaBegin(int nFiles, int nvars);
void printProcessedFiles(std::string fileName, int n);
void printProcessingFiles();
void printDone(int nbs);
void printOutputLocation1(std::string outFile);
void printOutputLocation2(std::string outFile);
void printOutputLocation3(std::string outFile);
void printZeroVariantsError(std::string fileName);
void printNColumnError(int z, int n, std::string fileName, int nheader);
void printOpenFileError(std::string fileName);
void printTimeCompleted(double wall0, double wall1, double cpu0, double cpu1);
#endif
