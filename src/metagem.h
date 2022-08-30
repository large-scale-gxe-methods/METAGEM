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

#define VERSION "1.0"

using std::cout;
using std::endl;
using std::cerr;
using spp::sparse_hash_map;

int flipAllele(std::string a1, std::string a2, std::string b1, std::string b2);
void metagem(CommandLine cmd);


#endif
