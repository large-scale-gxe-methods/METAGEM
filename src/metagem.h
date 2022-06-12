#ifndef METAGEM_H
#define METAGEM_H

#include <iostream>
#include <fstream>
#include <chrono>
#include <unordered_map>

#include <boost/math/distributions/chi_squared.hpp>

#include "../include/CLI11-2.2.0/CLI11.hpp"
#include "../include/sparsepp-2018.2/spp.h"
#include "matrix_utils.h"
#include "read_parameters.h"
#include "file.h"
#include "print.h"

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

#endif
