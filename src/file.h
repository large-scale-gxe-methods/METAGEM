#ifndef FILE_H
#define FILE_H

class FileInfo {
    public:

        int fileCount = 0;
        int nBetas = 0;

        std::vector<std::string> fileNames;
        std::vector<std::string> lc_betaIntNames;
        std::vector<std::string> betaIntNames;


        // Indices
        std::unordered_map<std::string, int> nvars;
        std::unordered_map<std::string, int> snpColumn;
        std::unordered_map<std::string, int> chrColumn;
        std::unordered_map<std::string, int> posColumn;
        std::unordered_map<std::string, int> effectColumn;
        std::unordered_map<std::string, int> nonEffectColumn;
        std::unordered_map<std::string, int> nSampleColumn;
        std::unordered_map<std::string, int> freqColumn;

        std::unordered_map<std::string, int> betaMargColumn;
        std::unordered_map<std::string, int> mb_seMargColumn;
        std::unordered_map<std::string, int> rb_seMargColumn;

        std::unordered_map<std::string, std::vector<int>> betaIntColumn;
        std::unordered_map<std::string, std::vector<int>> mb_covIntColumn;
        std::unordered_map<std::string, std::vector<int>> rb_covIntColumn;      
};

void processFileHeader(int Sq1, int metaOpt, std::vector<std::string> lc_intNames, std::string fileName, FileInfo* fip);

void printOutputHeader(std::string output, int metaOpt, size_t Sq1, std::vector<std::string> intNames); 


#endif
