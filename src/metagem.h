#ifndef METAGEM_H
#define METAGEM_H

#include <iostream>
#include <fstream>
#include <chrono>
#include <unordered_map>

#include <boost/math/distributions/chi_squared.hpp>
#include <boost/filesystem.hpp>

#include "matrix_utils.h"
#include "read_parameters.h"
#include "file.h"

using std::cout;
using std::endl;
using std::cerr;

int flipAllele(std::string a1, std::string a2, std::string b1, std::string b2);
void metagem(CommandLine cmd);



#endif