// <METAGEM: META-analysis of GEM summary statistics>
// Copyright (C) <2021-2023> Duy T. Pham and Han Chen 
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.


#include "metagem.h"
#include "time.h"

int main(int argc, char* argv[]) 
{
    printWelcome();
    
    // Start Timers
    double wall0 = get_wall_time();
    double cpu0  = get_cpu_time();

    // Process command line arguments
    CommandLine cmd;
    cmd.processCommandLine(argc, argv);

    metagem(cmd);

    // Stop timers
    double wall1 = get_wall_time();
    double cpu1  = get_cpu_time();
    printTimeCompleted(wall0, wall1, cpu0, cpu1);
    
    return(0);
}


void metagem(CommandLine cmd)
{
    // Read the header of each file
    printProcessingFiles();
    FileInfo* fip = new FileInfo();
    processFileHeader(cmd.nInt + 1, cmd.nInt2, cmd.mb, cmd.rb, cmd.lcIntNames, cmd.lcIntNames2, cmd.fileNames, fip);

    // Initialize objects
    int nvars = 0;
    bool mb = cmd.mb;
    bool rb = cmd.rb;
    size_t nInt   = cmd.nInt;
    size_t nInt1  = nInt + 1; // plus G
    size_t nFiles = cmd.fileNames.size();

    int nInt1_sq = nInt1 * nInt1;
    std::vector<std::string> snpid;
    std::vector<std::string> effectAllele;
    std::vector<std::string> nonEffectAllele;
    std::vector<double> nSamples;
    std::vector<double> AF;
    std::vector<size_t> nSeen;

    std::vector<double> mb_MU;
    std::vector<double> mb_MV;
    std::vector<double> mb_U;
    std::vector<double> mb_V;        
    std::vector<double> rb_MU;
    std::vector<double> rb_MV;
    std::vector<double> rb_U;
    std::vector<double> rb_V;
    sparse_hash_map<std::pair<std::string, std::string>, int> snpid_idx;
    
    for (size_t f = 0; f < nFiles; f++)
    {
        // Initialize
        std::string fileName = fip->fileNames[f];
        int snpColumn        = fip->snpColumn[fileName];
        int chrColumn        = fip->chrColumn[fileName];
        int posColumn        = fip->posColumn[fileName];
        int effectColumn     = fip->effectColumn[fileName];
        int nonEffectColumn  = fip->nonEffectColumn[fileName];
        int nSampleColumn    = fip->nSampleColumn[fileName];
        int afColumn         = fip->freqColumn[fileName];
        int betaMargColumn   = fip->betaMargColumn[fileName];	
        int nheader          = fip->nheader[fileName];
        std::vector<int> betaIntColumn = fip->betaIntColumn[fileName];

        // Read input file
        std::ifstream file;
        file.open(fileName);
        if (!file.is_open()) {
            printOpenFileError(fileName);
            exit(1);
        }

        // Ignore comments in the text file (ex. #dispersion: )
        std::string line;
        getline(file, line);

        // Read in summary statistics
        int n = 0;
        std::string value;
        std::vector <std::string> values(nheader);
        while(line.rfind("#", 0) == 0) { getline(file, line); }
    
        
        std::vector<int> seen(nvars, 0);
        while(getline(file, line))
        {
            line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
            
            // Get values
            int z = 0;
            std::istringstream iss(line);
            while (getline(iss, value, '\t')) {
                values[z] = value;
                z++;
            }

            if (z != nheader) {
                printNColumnError(z, n + 1, fileName, nheader);
                exit(1);
            }

            std::string snpName = values[snpColumn];
            std::string chr = values[chrColumn];
            std::string pos = values[posColumn];
            std::string a1  = values[nonEffectColumn];
            std::string a2  = values[effectColumn];
            double file_ns  = std::stod(values[nSampleColumn]);
            double file_af  = std::stod(values[afColumn]);
            std::transform(a1.begin(), a1.end(), a1.begin(), [](char c){ return std::tolower(c); });
            std::transform(a2.begin(), a2.end(), a2.begin(), [](char c){ return std::tolower(c); });

            
            std::string identifier = chr + ":" + pos;
            std::string idA = identifier + ":" + a1 + ":" + a2;
            std::string idB = identifier + ":" + a2 + ":" + a1;
            std::pair<std::string, std::string> id = {idA, idB};
            sparse_hash_map<std::pair<std::string, std::string>, int>::iterator it = snpid_idx.find(id);
            if (it == snpid_idx.end()) {
                snpid_idx[id] = nvars;
                seen.push_back(1);

                snpid.push_back(snpName + "\t" + chr+ "\t" + pos + "\t" + a1 + "\t" + a2 + "\t");
                nonEffectAllele.push_back(a1);
                effectAllele.push_back(a2);
                AF.push_back(file_af * 2.0 * file_ns);
                nSamples.push_back(file_ns);
                nSeen.push_back(1);
                
                // Matrix indices
                int ui = nvars * nInt1;
                int vi = nvars * nInt1_sq;

                // Betas
                double fileBM = std::stod(values[betaMargColumn]);
                std::vector<double> fileBetaInt(nInt1);
                for (size_t i = 0; i < nInt1; i++) {
                    fileBetaInt[i] = std::stod(values[betaIntColumn[i]]);
                }

                // Model-based
                if (mb)
                {
                    // Marginal
                    int mb_seMargColumn = fip->mb_seMargColumn[fileName];
                    double mb_fileMV    = std::stod(values[mb_seMargColumn]);
                    double mb_margV     = 1 / (mb_fileMV * mb_fileMV);
                    mb_MV.push_back(1 / (mb_fileMV * mb_fileMV));
                    mb_MU.push_back(mb_margV * (fileBM));

                    // Interaction
                    mb_U.resize(mb_U.size() + nInt1);
                    std::vector<int> mb_covIntColumn = fip->mb_covIntColumn[fileName];
                    for (size_t i = 0; i < nInt1; i++) {
                        size_t isq = (i * nInt1);
                        size_t ii  = vi + isq;
                        for (size_t j = 0; j < nInt1; j++) {
                            mb_V.push_back(std::stod(values[mb_covIntColumn[isq + j]]));
                        }
                        mb_V[ii + i] *= mb_V[ii + i];
                    }

                    // Add to main matrix
                    subMatInv(&mb_V[0], nInt1, vi);
                    subMatsubVecprod(&mb_V[0], &fileBetaInt[0], &mb_U[0], nInt1, nInt1, vi, 0, ui);
                }

                // Robust
                if (rb)
                {
                    // Marginal
                    int rb_seMargColumn = fip->rb_seMargColumn[fileName];
                    double rb_fileMV    = std::stod(values[rb_seMargColumn]);
                    double rb_margV     = 1 / (rb_fileMV * rb_fileMV);
                    rb_MV.push_back(1 / (rb_fileMV * rb_fileMV));
                    rb_MU.push_back(rb_margV * (fileBM));

                    // Interaction
                    rb_U.resize(rb_U.size() + nInt1);
                    std::vector<int> rb_covIntColumn = fip->rb_covIntColumn[fileName];
                    for (size_t i = 0; i < nInt1; i++) {
                        size_t isq = (i * nInt1);
                        size_t ii  = vi + isq;
                        for (size_t j = 0; j < nInt1; j++) {
                            rb_V.push_back(std::stod(values[rb_covIntColumn[isq + j]]));
                        }
                        rb_V[ii + i] *= rb_V[ii + i];
                    }

                    // Add to main matrix
                    subMatInv(&rb_V[0], nInt1, vi);
                    subMatsubVecprod(&rb_V[0], &fileBetaInt[0], &rb_U[0], nInt1, nInt1, vi, 0, ui);
                }

                nvars++;
            } else {
                
                int index = it->second;
                if (seen[index] == 0) {
                    seen[index] = 1;
                } else {
                    n++;
                    continue;
                }

                int dir = 1.0;
                std::string na = nonEffectAllele[index];
                std::string ea = effectAllele[index];
                if ((a1.compare(na) == 0) && (a2.compare(ea) == 0)) {
                    AF[index] += file_af * 2.0 * file_ns;

                } else if ((a1.compare(ea) == 0) && (a2.compare(na) == 0)){
                    dir = -1.0;
                    AF[index] += (1 - file_af) * 2.0 * file_ns;
                }

                nSeen[index] += 1;
                nSamples[index] += file_ns;


                // Matrix index
                int ui = index * nInt1;
                int vi = index * (nInt1 * nInt1);

                // Betas
                double fileBM = std::stod(values[betaMargColumn]) * dir;
                std::vector<double> fileBetaInt(nInt1);
                for (size_t i = 0; i < nInt1; i++) {
                    fileBetaInt[i] = std::stod(values[betaIntColumn[i]]) * dir;
                }

                // Model-based
                if (mb)
                {
                    // Marginal
                    int mb_seMargColumn = fip->mb_seMargColumn[fileName];
                    double mb_fileMV = std::stod(values[mb_seMargColumn]);

                    mb_fileMV = 1 / (mb_fileMV * mb_fileMV);
                    mb_MV[index] += mb_fileMV;
                    mb_MU[index] += (mb_fileMV * (fileBM));

                    // Interactions
                    std::vector<double> mb_fileU(nInt1);
                    std::vector<double> mb_fileV(nInt1 * nInt1);
                    std::vector<int> mb_covIntColumn = fip->mb_covIntColumn[fileName];

                    for (size_t i = 0; i < nInt1; i++) {
                        size_t ii = i * nInt1;
                        for (size_t j = 0; j < nInt1; j++) {
                            mb_fileV[ii + j] = std::stod(values[mb_covIntColumn[ii + j]]);
                        }
                        mb_fileV[ii + i] *= mb_fileV[ii + i];
                    }

                    // Add to main matrix
                    matInv(&mb_fileV[0], nInt1);
                    matvecprod(&mb_fileV[0], &fileBetaInt[0], &mb_fileU[0], nInt1, nInt1);
                    for (size_t i = 0; i < nInt1; i++){
                        size_t ii = i * nInt1;
                        mb_U[ui + i] += mb_fileU[i];
                        for (size_t j = 0; j < nInt1; j++) {
                            mb_V[vi + ii + j] += mb_fileV[ii + j];
                        }
                    }
                }

                // Robust
                if (rb)
                {
                    // Marginal
                    int rb_seMargColumn = fip->rb_seMargColumn[fileName];
                    double rb_fileMV = std::stod(values[rb_seMargColumn]);

                    rb_fileMV = 1 / (rb_fileMV * rb_fileMV);
                    rb_MV[index] += rb_fileMV;
                    rb_MU[index] += (rb_fileMV * (fileBM));

                    // Interactions
                    std::vector<double> rb_fileU(nInt1);
                    std::vector<double> rb_fileV(nInt1 * nInt1);
                    std::vector<int> rb_covIntColumn = fip->rb_covIntColumn[fileName];
                    
                    for (size_t i = 0; i < nInt1; i++) {
                        size_t ii = i * nInt1;
                        for (size_t j = 0; j < nInt1; j++) {
                            rb_fileV[ii + j] = std::stod(values[rb_covIntColumn[ii + j]]);
                        }
                        rb_fileV[ii + i] *= rb_fileV[ii + i];
                    }

                    // Add to main matrix
                    matInv(&rb_fileV[0], nInt1);
                    matvecprod(&rb_fileV[0], &fileBetaInt[0], &rb_fileU[0], nInt1, nInt1);
                    for (size_t i = 0; i < nInt1; i++) {
                        size_t ii = i * nInt1;
                        rb_U[ui + i] += rb_fileU[i];
                        for (size_t j = 0; j < nInt1; j++) {
                            rb_V[vi + ii + j] += rb_fileV[ii + j];
                        }
                    }
                }
            }
            n++;
        }
        file.close();

        if (n == 0) {
            printZeroVariantsError(fileName);
            exit(1);
        }

        printProcessedFiles(fileName, n);
    }
    printDone(2);

    // Create the output file and write column header names
    printOutputHeader(mb, rb, cmd.outFile, nInt1, cmd.intNames);
    std::ofstream results(cmd.outFile, std::ios_base::app);
    std::ostringstream oss;




    double vMarg;
    double betaMarg;
    double varMarg;
    double statMarg;
    double statInt;
    double statJoint;
    double pvalMarg;
    double pvalInt;
    double pvalJoint;
    double* Ai = new double[nInt1 * nInt1];
    double* VE = new double[nInt * nInt];
    std::vector<double> StempE(nInt, 0.0);
    std::vector<double> StempGE(nInt1, 0.0);
    std::vector<double> betaInt(nInt1, 0.0);
    boost::math::chi_squared chisq_dist_M(1);
    boost::math::chi_squared chisq_dist_Int(nInt);
    boost::math::chi_squared chisq_dist_Joint(nInt1);

    printMetaBegin(nFiles, nvars);
    for (int i = 0; i < nvars; i++)
    {
        int is  = (i * nInt1);
        int iss = (i * nInt1 * nInt1);

        oss << snpid[i] << nSeen[i] << "\t" << nSamples[i] << "\t" << AF[i] / 2.0 / nSamples[i] << "\t";
        if (mb)
        {
            subMatrix(&mb_V[0], Ai, nInt1, nInt1, nInt1, nInt1, iss);
            subMatInv(&mb_V[0], nInt1, iss);
            subMatrix(&mb_V[0], VE, nInt, nInt, nInt1, nInt, iss + nInt1 + 1);

            // Model-based Marginal Variance
            vMarg = mb_MV[i];
            varMarg  = 1 / vMarg;
            betaMarg = mb_MU[i] * varMarg;
            statMarg = (betaMarg * betaMarg) * vMarg;
            pvalMarg = (std::isnan(statMarg) || statMarg <= 0.0) ? NAN : boost::math::cdf(complement(chisq_dist_M, statMarg));


            // Interaction effects
            for (size_t j = 0; j < nInt1; j++) {
                for (size_t k = 0; k < nInt1; k++) {
                    betaInt[j] += (mb_V[iss + (nInt1 * j) + k] * mb_U[is + k]);
                }
            }


            // Int P-value
            matInv(VE, nInt);
            for (size_t j = 0; j < nInt; j++) {
                for (size_t k = 0; k < nInt; k++) {
                    StempE[j] += (VE[(nInt * j) + k] * betaInt[k + 1]);
                }
            }
            
            statInt = 0.0;
            for (size_t j = 1; j < nInt1; j++) 
                statInt += betaInt[j] * StempE[j-1];
            pvalInt = (std::isnan(statInt) || statInt <= 0.0) ? NAN : boost::math::cdf(complement(chisq_dist_Int, statInt));
    

            // Joint P-value
            for (size_t j = 0; j < nInt1; j++) {
                for (size_t k = 0; k < nInt1; k++) {
                    StempGE[j] += (Ai[(nInt1 * j) + k] * betaInt[k]);
                }
            }

            statJoint = 0.0;
            for (size_t k = 0; k < nInt1; k++)
                statJoint += betaInt[k] * StempGE[k];
            pvalJoint = (std::isnan(statJoint) || statJoint <= 0.0) ? NAN : boost::math::cdf(complement(chisq_dist_Joint, statJoint));


            // Print
            oss << betaMarg << "\t" << sqrt(varMarg) << "\t";
            for (size_t j = 0; j < nInt1; j++) {
                oss << betaInt[j] << "\t";
            }
            for (size_t ii = 0; ii < nInt1; ii++) {
                oss << sqrt(mb_V[iss + (ii * nInt1) + ii]) << "\t";
            }
            for (size_t ii = 0; ii < nInt1; ii++) {
                for (size_t jj = 0; jj < nInt1; jj++) {
                    if (ii < jj) {
                        oss << mb_V[iss + (ii * nInt1) + jj] << "\t";
                    }
                }
            }
            oss << pvalMarg << "\t" << pvalInt << "\t" << pvalJoint << ((rb) ? "\t" : "\n");

            std::fill(StempE.begin(), StempE.end(), 0.0);
            std::fill(StempGE.begin(), StempGE.end(), 0.0);
            std::fill(betaInt.begin(), betaInt.end(), 0.0);
        }

        if (rb)
        {
            subMatrix(&rb_V[0], Ai, nInt1, nInt1, nInt1, nInt1, iss);
            subMatInv(&rb_V[0], nInt1, iss);
            subMatrix(&rb_V[0], VE, nInt, nInt, nInt1, nInt, iss + nInt1 + 1);

            // Robust Marginal Variance
            vMarg = rb_MV[i];
            varMarg  = 1 / vMarg;
            betaMarg = rb_MU[i] * varMarg;
            statMarg = (betaMarg * betaMarg) * vMarg;
            pvalMarg = (std::isnan(statMarg) || statMarg <= 0.0) ? NAN : boost::math::cdf(complement(chisq_dist_M, statMarg));


            // Interaction effects
            for (size_t j = 0; j < nInt1; j++) {
                for (size_t k = 0; k < nInt1; k++) {
                    betaInt[j] += (rb_V[iss + (nInt1 * j) + k] * rb_U[is + k]);
                }
            }


            // Int P-value
            matInv(VE, nInt);
            for (size_t j = 0; j < nInt; j++) {
                for (size_t k = 0; k < nInt; k++) {
                    StempE[j] += (VE[(nInt * j) + k] * betaInt[k + 1]);
                }
            }

            statInt = 0.0;
            for (size_t j = 1; j < nInt1; j++) 
                statInt += betaInt[j] * StempE[j-1];
            pvalInt = (std::isnan(statInt) || statInt <= 0.0) ? NAN : boost::math::cdf(complement(chisq_dist_Int, statInt));
    

            // Joint P-value
            for (size_t j = 0; j < nInt1; j++) {
                for (size_t k = 0; k < nInt1; k++) {
                    StempGE[j] += (Ai[(nInt1 * j) + k] * betaInt[k]);
                }
            }

            statJoint = 0.0;
            for (size_t k = 0; k < nInt1; k++)
                statJoint += betaInt[k] * StempGE[k];
            pvalJoint = (std::isnan(statJoint) || statJoint <= 0.0) ? NAN : boost::math::cdf(complement(chisq_dist_Joint, statJoint));


            // Print
            oss << betaMarg << "\t" << sqrt(varMarg) << "\t";
            for (size_t j = 0; j < nInt1; j++) {
                oss << betaInt[j] << "\t";
            }
            for (size_t ii = 0; ii < nInt1; ii++) {
                oss << sqrt(rb_V[iss + (ii * nInt1) + ii]) << "\t";
            }
            for (size_t ii = 0; ii < nInt1; ii++) {
                for (size_t jj = 0; jj < nInt1; jj++) {
                    if (ii < jj) {
                        oss << rb_V[iss + (ii * nInt1) + jj] << "\t";
                    }
                }
            }
            oss << pvalMarg << "\t" << pvalInt << "\t" << pvalJoint << "\n";

            std::fill(StempE.begin(), StempE.end(), 0.0);
            std::fill(StempGE.begin(), StempGE.end(), 0.0);
            std::fill(betaInt.begin(), betaInt.end(), 0.0);
        }

        if (i % 100000 == 0)
        {
            results << oss.str();
            oss.str(std::string());
            oss.clear();
        }
    }

    results << oss.str();
    oss.str(std::string());
    oss.clear();
    results.close();

    delete[] Ai;
    delete[] VE;

    printDone(2);
    printOutputLocation(cmd.outFile);
}

void printWelcome() {
    cout << "\n*********************************************************\n";
    cout << "Welcome to METAGEM v" << VERSION << "\n";
    cout << "(C) 2021-2023 Duy Pham and Han Chen \n";
    cout << "GNU General Public License v3\n";
    cout << "*********************************************************\n";
}

void printMetaBegin(int nFiles, int nvars) {
  cout << "Performing meta-analysis with " << std::to_string(nFiles) << " files and " << std::to_string(nvars) << " variants...\n";
}

void printProcessedFiles(std::string fileName, int n) {
  cout << "Processed file [" << fileName << "] with " << n << " variants.\n";
}

void printProcessingFiles() {
  cout << "Processing files...\n";
}

void printDone(int nbs) {

  if (nbs == 1) {
    cout << "Done.\n";
  } else if (nbs == 2) {
    cout << "Done.\n\n";
  } else {
    cerr << "\n ERROR: Invalid number.\n\n";
  }
}

void printOutputLocation(std::string outFile) {
  cout << "Results are in [" << outFile << "].";
}

void printZeroVariantsError(std::string fileName) {
  cerr << "\nERROR: No variants in file [" << fileName << "].\n\n";
}

void printNColumnError(int z, int n, std::string fileName, int nheader) {
  cerr << "\nERROR: Unexpect number of columns (" << z << ") for variant number " << n+1 << "in file [" << fileName << "]. Expected: " << nheader << ".\n\n";
}

void printOpenFileError(std::string fileName) {
    cerr << "\nERROR: Cannot open the file [" << fileName << "].\n\n";
}

void printTimeCompleted(double wall0, double wall1, double cpu0, double cpu1) {
    cout << "\n*********************************************************\n";
    cout << "Wall Time = " << wall1 - wall0 << " (sec)\n";
    cout << "CPU Time  = " << cpu1  - cpu0  << " (sec)\n\n";
}
