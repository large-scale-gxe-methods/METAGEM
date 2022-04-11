#include "metagem.h"
#include "time.h"

int main(int argc, char* argv[]) 
{

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
    cout << "\n*********************************************************\n";
    cout << "Wall Time = " << wall1 - wall0 << " (sec)\n";
    cout << "CPU Time  = " << cpu1  - cpu0  << " (sec)\n\n";
    
    return(0);
}


int flipAllele(std::string a1, std::string a2, std::string b1, std::string b2)
{
    if ((a1.compare(b1) == 0) && (a2.compare(b2) == 0)) {
        return 0;

    }
    else if ((a1.compare(b2) == 0) && (a2.compare(b1) == 0)){
        return 1;

    } else {
        return 2;
    }
}



void metagem(CommandLine cmd)
{
    bool mb = cmd.mb;
    bool rb = cmd.rb;

    int Sq1 = 1;
    size_t nFiles = cmd.fileNames.size();
    size_t nExp   = cmd.nExp;
    size_t nExp1  = nExp + 1; // plus G

    std::vector<std::string> lc_intNames = cmd.intNames;
    for(std::string &s : lc_intNames){
        std::transform(s.begin(), s.end(), s.begin(), [](char c){ return std::tolower(c); });
        s = "g-" + s;
        Sq1++;
    }
    lc_intNames.insert(lc_intNames.begin(), "g");


    // Read the header of each file
    FileInfo* fip = new FileInfo();
    processFileHeader(Sq1, cmd.metaOpt, lc_intNames, cmd.fileNames, fip);
    

    // Store unique variant indices
    int nvars = 0;
    sparse_hash_map<std::string, int> snpid_idx;

    std::vector<std::string> snpid;
    std::vector<std::string> effectAllele;
    std::vector<std::string> nonEffectAllele;
    std::vector<double> nSamples;
    std::vector<double> AF;
    std::vector<size_t> nfiles;

    std::vector<double> mb_MU;
    std::vector<double> rb_MU;
    std::vector<double> mb_MV;
    std::vector<double> rb_MV;
    std::vector<double> mb_U;
    std::vector<double> rb_U;
    std::vector<double> mb_V;
    std::vector<double> rb_V;

    boost::math::chi_squared chisq_dist_M(1);
    boost::math::chi_squared chisq_dist_Int(nExp);
    boost::math::chi_squared chisq_dist_Joint(nExp1);

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
        std::vector<int> betaIntColumn = fip->betaIntColumn[fileName];

        sparse_hash_map<std::string, int> seen;

        // Read input file
        std::ifstream file;
        file.open(fileName);
        if (!file.is_open()) {
            cerr << "\nERROR: Cannot open the file [" << fileName << "].\n\n";
            exit(1);
        }

        // Ignore comments in the text file (ex. #dispersion: )
        std::string line;
        getline(file, line);
        while(line.rfind("#", 0) == 0) { getline(file, line); }


        // Read in summary statistics
        int n = 0;
        while(getline(file, line))
        {
            line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
            
            // Get values
            std::istringstream iss(line);
            std::string value;
            std::vector <std::string> values;
            while (getline(iss, value, '\t')) values.push_back(value);


            std::string snpName = values[snpColumn];

            std::string a1 = values[nonEffectColumn];
            std::string a2 = values[effectColumn];
            std::transform(a1.begin(), a1.end(), a1.begin(), [](char c){ return std::tolower(c); });
            std::transform(a2.begin(), a2.end(), a2.begin(), [](char c){ return std::tolower(c); });

            double file_ns = std::stod(values[nSampleColumn]);


            sparse_hash_map<std::string, int>::iterator it = snpid_idx.find(snpName);
            if (it == snpid_idx.end()) {
                snpid_idx[snpName] = nvars;
                seen[snpName] = 1;

                snpid.push_back(snpName + "\t" + values[chrColumn] + "\t" + values[posColumn] + "\t" + a1 + "\t" + a2 + "\t");
                nonEffectAllele.push_back(a1);
                effectAllele.push_back(a2);
                AF.push_back(std::stod(values[afColumn]) * 2.0 * file_ns);
                nSamples.push_back(file_ns);
                nfiles.push_back(1);
                
                // Matrix indices
                int vi = nvars * (Sq1 * Sq1);

                // Betas
                double fileBM = std::stod(values[betaMargColumn]);
                std::vector<double> fileBetaInt(Sq1);
                for (int i = 0; i < Sq1; i++) {
                    fileBetaInt[i] = std::stod(values[betaIntColumn[i]]);
                }

                // Model-based
                if (mb)
                {
                    // Marginal
                    int mb_seMargColumn = fip->mb_seMargColumn[fileName];
                    double mb_fileMV    = std::stod(values[mb_seMargColumn]);

                    mb_MV.push_back(1 / (mb_fileMV * mb_fileMV));
                    mb_MU.push_back(mb_MV[nvars] * (fileBM));

                    // Interaction
                    std::vector<int> mb_covIntColumn = fip->mb_covIntColumn[fileName];
                    for (int i = 0; i < Sq1; i++) {
                        int isq = (i * Sq1);
                        int ii  = vi + isq;
                        for (int j = 0; j < Sq1; j++) {
                            mb_V.push_back(std::stod(values[mb_covIntColumn[isq + j]]));
                        }
                        mb_V[ii + i] *= mb_V[ii + i];
                    }

                    // Add to main matrix
                    std::vector<double> mb_fileU(Sq1);
                    subMatInv(&mb_V[0], Sq1, vi);
                    subMatsubVecprod(&mb_V[0], &fileBetaInt[0], &mb_fileU[0], Sq1, Sq1, vi, 0, 0);
                    
                    for (int i = 0; i < Sq1; i++) {
                        mb_U.push_back(mb_fileU[i]);
                    }
                }

                // Robust
                if (rb)
                {
                    // Marginal
                    int rb_seMargColumn = fip->rb_seMargColumn[fileName];
                    double rb_fileMV    = std::stod(values[rb_seMargColumn]);

                    rb_MV.push_back(1 / (rb_fileMV * rb_fileMV));
                    rb_MU.push_back(rb_MV[nvars] * (fileBM));

                    // Interaction
                    std::vector<int> rb_covIntColumn = fip->rb_covIntColumn[fileName];
                    for (int i = 0; i < Sq1; i++) {
                        int isq = (i * Sq1);
                        int ii  = vi + isq;
                        for (int j = 0; j < Sq1; j++) {
                            rb_V.push_back(std::stod(values[rb_covIntColumn[isq + j]]));
                        }
                        rb_V[ii + i] *= rb_V[ii + i];
                    }

                    // Add to main matrix
                    std::vector<double> rb_fileU(Sq1);
                    subMatInv(&rb_V[0], Sq1, vi);
                    subMatsubVecprod(&rb_V[0], &fileBetaInt[0], &rb_fileU[0], Sq1, Sq1, vi, 0, 0);

                    for (int i = 0; i < Sq1; i++) {
                        rb_U.push_back(rb_fileU[i]);
                    }
                }

                nvars++;
            } else {
                sparse_hash_map<std::string, int>::iterator it2 = seen.find(snpName);
                if (it2 == seen.end()) {
                    seen[snpName] = 1;
                } else {
                    n++;
                    continue;
                }

                int index = snpid_idx[snpName];

                nSamples[index] += file_ns;
                nfiles[index]++;

                // Matrix index
                int ui = index * Sq1;
                int vi = index * (Sq1 * Sq1);

                int dir  = 1.0;
                int flip = flipAllele(nonEffectAllele[index], effectAllele[index], a1, a2);
                double file_af = std::stod(values[afColumn]);
                if (flip == 0) {
                    AF[index] += file_af * 2.0 * file_ns;

                } else if (flip == 1) {
                    dir = -1.0;
                    AF[index] += (1 - file_af) * 2.0 * file_ns;

                } else {
                    n++;
                    continue;
                }

                // Betas
                double fileBM = std::stod(values[betaMargColumn]) * dir;
                std::vector<double> fileBetaInt(Sq1);
                for (int i = 0; i < Sq1; i++) {
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
                    std::vector<double> mb_fileU(Sq1);
                    std::vector<double> mb_fileV(Sq1 * Sq1);
                    std::vector<int> mb_covIntColumn = fip->mb_covIntColumn[fileName];

                    for (int i = 0; i < Sq1; i++) {
                        int ii = i * Sq1;
                        for (int j = 0; j < Sq1; j++) {
                            mb_fileV[ii + j] = std::stod(values[mb_covIntColumn[ii + j]]);
                        }
                        mb_fileV[ii + i] *= mb_fileV[ii + i];
                    }

                    // Add to main matrix
                    matInv(&mb_fileV[0], Sq1);
                    matvecprod(&mb_fileV[0], &fileBetaInt[0], &mb_fileU[0], Sq1, Sq1);
                    for (int i = 0; i < Sq1; i++){
                        int ii = i * Sq1;
                        mb_U[ui + i] += mb_fileU[i];
                        for (int j = 0; j < Sq1; j++) {
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
                    std::vector<double> rb_fileU(Sq1);
                    std::vector<double> rb_fileV(Sq1 * Sq1);
                    std::vector<int> rb_covIntColumn = fip->rb_covIntColumn[fileName];

                    for (int i = 0; i < Sq1; i++) {
                        int ii = i * Sq1;
                        for (int j = 0; j < Sq1; j++) {
                            rb_fileV[ii + j] = std::stod(values[rb_covIntColumn[ii + j]]);
                        }
                        rb_fileV[ii + i] *= rb_fileV[ii + i];
                    }

                    // Add to main matrix
                    matInv(&rb_fileV[0], Sq1);
                    matvecprod(&rb_fileV[0], &fileBetaInt[0], &rb_fileU[0], Sq1, Sq1);
                    for (int i = 0; i < Sq1; i++) {
                        int ii = i * Sq1;
                        rb_U[ui + i] += rb_fileU[i];
                        for (int j = 0; j < Sq1; j++) {
                            rb_V[vi + ii + j] += rb_fileV[ii + j];
                        }
                    }
                }
            } 
            n++;
        }
        file.close();

        if (n == 0) {
            cerr << "\nERROR: No variants in file [" << fileName << "].\n\n";
            exit(1);
        }

        cout << "Processed file [" << fileName << "] with " << std::to_string(n) << " variants.\n";
    }

    // Create the output file and write column header names
    printOutputHeader(cmd.outFile, cmd.metaOpt, Sq1, cmd.intNames);
    std::ofstream results(cmd.outFile, std::ios_base::app);
    std::ostringstream oss;


    double mb_betaMarg;
    double mb_varMarg;
    double rb_betaMarg;
    double rb_varMarg;
    double mb_pvalMarg;
    double mb_pvalInt;
    double mb_pvalJoint;
    double rb_pvalMarg;
    double rb_pvalInt;
    double rb_pvalJoint;

    double* mb_betaInt = new double[Sq1];
    double* rb_betaInt = new double[Sq1];
    double* StempE  = new double[nExp];
    double* StempGE = new double[nExp1];
    double* Ai      = new double[nExp1 * nExp1];
    double* VE      = new double[nExp * nExp];

    cout << "\nPerforming meta-analysis with " << std::to_string(nFiles) << " files and " << std::to_string(nvars) << " variants...\n";
    for (int i = 0; i < nvars; i++)
    {

        int is  = (i * Sq1);
        int iss = (i * Sq1 * Sq1);

        AF[i] = AF[i] / 2.0 / nSamples[i];

        if (mb)
        {
            subMatInv(&mb_V[0], Sq1, iss);

            // Model-based Marginal Variance
            mb_betaMarg = mb_MU[i] / mb_MV[i];
            mb_varMarg  = 1 / mb_MV[i];
            double mb_statMarg = (mb_betaMarg * mb_betaMarg) / mb_varMarg;
            mb_pvalMarg = (std::isnan(mb_statMarg) || mb_statMarg <= 0.0) ? NAN : boost::math::cdf(complement(chisq_dist_M, mb_statMarg));

            // Interaction effects
            subMatsubVecprod(&mb_V[0], &mb_U[0], mb_betaInt, Sq1, Sq1, iss, is, 0);

            // Model-based Int P-value
            subMatrix(&mb_V[0], VE, nExp, nExp, Sq1, nExp, iss + Sq1 + 1);
            matInv(VE, 1);
            matvecSprod(VE, mb_betaInt, StempE, nExp, nExp, 1);
            double mb_statInt = 0.0;
            for (size_t j = 1; j < nExp1; j++) 
                mb_statInt += mb_betaInt[j] * StempE[j-1];
            mb_pvalInt = (std::isnan(mb_statInt) || mb_statInt <= 0.0) ? NAN : boost::math::cdf(complement(chisq_dist_Int, mb_statInt));
    
            // Model-based Joint P-value
            subMatrix(&mb_V[0], Ai, nExp1, nExp1, nExp1, nExp1, iss);
            matInv(Ai, nExp1);
            matvecprod(Ai, mb_betaInt, StempGE, nExp1, nExp1);
            double mb_statJoint = 0.0;
            for (size_t k = 0; k < nExp1; k++)
                mb_statJoint += mb_betaInt[k] * StempGE[k];
            mb_pvalJoint = (std::isnan(mb_statJoint) || mb_statJoint <= 0.0) ? NAN : boost::math::cdf(complement(chisq_dist_Joint, mb_statJoint));
    
        }

        if (rb)
        {
            subMatInv(&rb_V[0], Sq1, iss);

            // Robust Marginal Variance
            rb_betaMarg = rb_MU[i] / rb_MV[i];
            rb_varMarg  = 1 / rb_MV[i];
            double rb_statMarg = (rb_betaMarg * rb_betaMarg) / rb_varMarg;
            rb_pvalMarg = (std::isnan(rb_statMarg) || rb_statMarg <= 0.0) ? NAN : boost::math::cdf(complement(chisq_dist_M, rb_statMarg));

            // Interaction effects
            subMatsubVecprod(&rb_V[0], &rb_U[0], rb_betaInt, Sq1, Sq1, iss, is, 0);

            // Robust Int P-value
            subMatrix(&rb_V[0], VE, nExp, nExp, Sq1, nExp, iss + Sq1 + 1);
            matInv(VE, 1);
            matvecSprod(VE, rb_betaInt, StempE, nExp, nExp, 1);
            double rb_statInt = 0.0;
            for (size_t j = 1; j < nExp1; j++) 
                rb_statInt += rb_betaInt[j] * StempE[j-1];
            rb_pvalInt = (std::isnan(rb_statInt) || rb_statInt <= 0.0) ? NAN : boost::math::cdf(complement(chisq_dist_Int, rb_statInt));

            // Robust Joint P-value
            subMatrix(&rb_V[0], Ai, nExp1, nExp1, nExp1, nExp1, iss);
            matInv(Ai, nExp1);
            matvecprod(Ai, rb_betaInt, StempGE, nExp1, nExp1);
            double rb_statJoint = 0.0;
            for (size_t k = 0; k < nExp1; k++)
                rb_statJoint += rb_betaInt[k] * StempGE[k];
            rb_pvalJoint = (std::isnan(rb_statJoint) || rb_statJoint <= 0.0) ? NAN : boost::math::cdf(complement(chisq_dist_Joint, rb_statJoint));
        }


        oss << snpid[i] << nfiles[i] << "\t" << nSamples[i] << "\t" << AF[i] << "\t";

        // Marginal summary stats
        if (mb)
        {
            oss << mb_betaMarg << "\t" << mb_varMarg << "\t";
        }
        if (rb)
        {
            oss << rb_betaMarg << "\t" << rb_varMarg << "\t";
        }

        // Model based summary stats
        if (mb) {
            for (int j = 0; j < Sq1; j++) {
                oss << mb_betaInt[j] << "\t";
            }

            for (int ii = 0; ii < Sq1; ii++) {
                oss << sqrt(mb_V[iss + (ii * Sq1) + ii]) << "\t";
            }

            for (int ii = 0; ii < Sq1; ii++) {
                for (int jj = 0; jj < Sq1; jj++) {
                    if (ii < jj) {
                        oss << mb_V[iss + (ii * Sq1) + jj] << "\t";
                    }
                }
            }
        }

        // Robust summary stats
        if (rb) 
        {
            for (int j = 0; j < Sq1; j++) {
                oss << rb_betaInt[j] << "\t";
            }

            for (int ii = 0; ii < Sq1; ii++) {
                oss << sqrt(rb_V[iss + (ii * Sq1) + ii]) << "\t";
            }

            for (int ii = 0; ii < Sq1; ii++) {
                for (int jj = 0; jj < Sq1; jj++) {
                    if (ii < jj) {
                        oss << rb_V[iss + (ii * Sq1) + jj] << "\t";
                    }
                }
            }
        }

        if (mb)
        {
            oss << mb_pvalMarg << "\t" << mb_pvalInt << "\t" << mb_pvalJoint << ((rb) ? "\t" : "\n");
        }

        if (rb)
        {
            oss << rb_pvalMarg << "\t" << rb_pvalInt << "\t" << rb_pvalJoint << "\n";
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

    delete[] mb_betaInt;
    delete[] rb_betaInt;
    delete[] StempE;
    delete[] StempGE;
    delete[] Ai;
    delete[] VE;

    cout << "Done.\n\n";
    cout << "Results are in [" << cmd.outFile << "].";
}