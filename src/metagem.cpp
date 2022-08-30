#include "metagem.h"
#include "time.h"
#include "print.h"

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
    processFileHeader(cmd.nInt + 1, cmd.mb, cmd.rb, cmd.lcIntNames, cmd.fileNames, fip);

    // Initialize objects
    int nvars = 0;
    bool mb = cmd.mb;
    bool rb = cmd.rb;
    size_t nInt   = cmd.nInt;
    size_t nInt1  = nInt + 1; // plus G
    size_t nFiles = cmd.fileNames.size();

    int init = fip->file0_nvars;
    int nInt1_sq = nInt1 * nInt1;
    std::vector<std::string> snpid(init);
    std::vector<std::string> effectAllele(init);
    std::vector<std::string> nonEffectAllele(init);
    std::vector<double> nSamples(init, 0.0);
    std::vector<double> AF(init, 0.0);
    std::vector<size_t> nSeen(init, 0);

    std::vector<double> mb_MU;
    std::vector<double> mb_MV;
    std::vector<double> mb_U;
    std::vector<double> mb_V;        
    std::vector<double> rb_MU;
    std::vector<double> rb_MV;
    std::vector<double> rb_U;
    std::vector<double> rb_V;
    
    if (mb) {
        mb_MU.resize(init);
        mb_MV.resize(init);
        mb_U.resize(init * nInt1);
        mb_V.resize(init * nInt1_sq);        
    }
    if (rb) {
        rb_MU.resize(init);
        rb_MV.resize(init);
        rb_U.resize(init * nInt1);
        rb_V.resize(init * nInt1_sq);
    }

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

        std::vector<uint8_t> seen;

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
    
        if (f == 0) {
            sparse_hash_map<std::string, std::vector<int>> seen;
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

                sparse_hash_map<std::string, std::vector<int>>::iterator it = seen.find(identifier);
                if (it != seen.end()) {
                    bool duplicate = false;
                    std::vector<int> indices = it->second;
                    for (size_t i = 0; i < indices.size(); i++) {
                        std::string na = nonEffectAllele[indices[i]];
                        std::string ea = effectAllele[indices[i]];
                        if (((a1.compare(na) == 0) && (a2.compare(ea) == 0)) ||  ((a1.compare(ea) == 0) && (a2.compare(na) == 0))) {
                            duplicate = true;
                        }
                    }

                    if (duplicate) {
                        n++;
                        continue;
                    }
                }
                seen[identifier].push_back(1);

                int index = nvars;
                snpid_idx[id] = index;
                snpid[index] = snpName + "\t" + chr+ "\t" + pos + "\t" + a1 + "\t" + a2 + "\t";
                nonEffectAllele[index] = a1;
                effectAllele[index]    = a2;
                AF[index] = file_af * 2.0 * file_ns;
                nSamples[index] = file_ns;
                nSeen[index] = 1;
                
                // Matrix index
                int ui = index * nInt1;
                int vi = index * (nInt1 * nInt1);

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

                nvars++;
                n++;
            }

            if (n != nvars) {
                snpid.resize(nvars);
                effectAllele.resize(nvars);
                nonEffectAllele.resize(nvars);
                nSamples.resize(nvars);
                AF.resize(nvars);
                nSeen.resize(nvars);

                if (mb) {
                    mb_MU.resize(nvars);
                    mb_MV.resize(nvars);
                    mb_U.resize(nvars * nInt1);
                    mb_V.resize(nvars * nInt1_sq);        
                }
                if (rb) {
                    rb_MU.resize(nvars);
                    rb_MV.resize(nvars);
                    rb_U.resize(nvars * nInt1);
                    rb_V.resize(nvars * nInt1_sq);
                }
            }
        } else {
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
    double* Ai = new double[nInt1 * nInt1];
    double* VE = new double[nInt * nInt];
    double* StempE  = new double[nInt];
    double* StempGE = new double[nInt1];
    double* mb_betaInt = new double[nInt1];
    double* rb_betaInt = new double[nInt1];

    boost::math::chi_squared chisq_dist_M(1);
    boost::math::chi_squared chisq_dist_Int(nInt);
    boost::math::chi_squared chisq_dist_Joint(nInt1);


    printMetaBegin(nFiles, nvars);
    for (int i = 0; i < nvars; i++)
    {

        int is  = (i * nInt1);
        int iss = (i * nInt1 * nInt1);

        AF[i] = AF[i] / 2.0 / nSamples[i];
        oss << snpid[i] << nSeen[i] << "\t" << nSamples[i] << "\t" << AF[i] << "\t";

        if (mb)
        {
            subMatrix(&mb_V[0], Ai, nInt1, nInt1, nInt1, nInt1, iss);
            subMatInv(&mb_V[0], nInt1, iss);

            // Model-based Marginal Variance
            double mb_vMarg = mb_MV[i];
            mb_varMarg  = 1 / mb_vMarg;
            mb_betaMarg = mb_MU[i] * mb_varMarg;
            double mb_statMarg = (mb_betaMarg * mb_betaMarg) * mb_vMarg;
            mb_pvalMarg = (std::isnan(mb_statMarg) || mb_statMarg <= 0.0) ? NAN : boost::math::cdf(complement(chisq_dist_M, mb_statMarg));

            // Interaction effects
            subMatsubVecprod(&mb_V[0], &mb_U[0], mb_betaInt, nInt1, nInt1, iss, is, 0);

            // Model-based Int P-value
            subMatrix(&mb_V[0], VE, nInt, nInt, nInt1, nInt, iss + nInt1 + 1);
            matInv(VE, nInt);
            matvecSprod(VE, mb_betaInt, StempE, nInt, nInt, 1);
            double mb_statInt = 0.0;
            for (size_t j = 1; j < nInt1; j++) 
                mb_statInt += mb_betaInt[j] * StempE[j-1];
            mb_pvalInt = (std::isnan(mb_statInt) || mb_statInt <= 0.0) ? NAN : boost::math::cdf(complement(chisq_dist_Int, mb_statInt));
    
            // Model-based Joint P-value
            matvecprod(Ai, mb_betaInt, StempGE, nInt1, nInt1);
            double mb_statJoint = 0.0;
            for (size_t k = 0; k < nInt1; k++)
                mb_statJoint += mb_betaInt[k] * StempGE[k];
            mb_pvalJoint = (std::isnan(mb_statJoint) || mb_statJoint <= 0.0) ? NAN : boost::math::cdf(complement(chisq_dist_Joint, mb_statJoint));

            // Print
            oss << mb_betaMarg << "\t" << sqrt(mb_varMarg) << "\t";
            for (size_t j = 0; j < nInt1; j++) {
                oss << mb_betaInt[j] << "\t";
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
            oss << mb_pvalMarg << "\t" << mb_pvalInt << "\t" << mb_pvalJoint << ((rb) ? "\t" : "\n");
    
        }

        if (rb)
        {
            subMatrix(&rb_V[0], Ai, nInt1, nInt1, nInt1, nInt1, iss);
            subMatInv(&rb_V[0], nInt1, iss);

            // Robust Marginal Variance
            double rb_vMarg = rb_MV[i];
            rb_varMarg  = 1 / rb_vMarg;
            rb_betaMarg = rb_MU[i] * rb_varMarg;
            double rb_statMarg = (rb_betaMarg * rb_betaMarg) * rb_vMarg;
            rb_pvalMarg = (std::isnan(rb_statMarg) || rb_statMarg <= 0.0) ? NAN : boost::math::cdf(complement(chisq_dist_M, rb_statMarg));

            // Interaction effects
            subMatsubVecprod(&rb_V[0], &rb_U[0], rb_betaInt, nInt1, nInt1, iss, is, 0);

            // Robust Int P-value
            subMatrix(&rb_V[0], VE, nInt, nInt, nInt1, nInt, iss + nInt1 + 1);
            matInv(VE, nInt);
            matvecSprod(VE, rb_betaInt, StempE, nInt, nInt, 1);
            double rb_statInt = 0.0;
            for (size_t j = 1; j < nInt1; j++) 
                rb_statInt += rb_betaInt[j] * StempE[j-1];
            rb_pvalInt = (std::isnan(rb_statInt) || rb_statInt <= 0.0) ? NAN : boost::math::cdf(complement(chisq_dist_Int, rb_statInt));

            // Robust Joint P-value
            matvecprod(Ai, rb_betaInt, StempGE, nInt1, nInt1);
            double rb_statJoint = 0.0;
            for (size_t k = 0; k < nInt1; k++)
                rb_statJoint += rb_betaInt[k] * StempGE[k];
            rb_pvalJoint = (std::isnan(rb_statJoint) || rb_statJoint <= 0.0) ? NAN : boost::math::cdf(complement(chisq_dist_Joint, rb_statJoint));

            // Print
            oss << rb_betaMarg << "\t" << sqrt(rb_varMarg) << "\t";
            for (size_t j = 0; j < nInt1; j++) {
                oss << rb_betaInt[j] << "\t";
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

    printDone(2);
    printOutputLocation(cmd.outFile);
}