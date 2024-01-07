#ifndef FILE_H
#define FILE_H

class FileInfo {
    public:

        int v = 0;
        int nBetas = 0;
        int fileCount = 0;

        std::vector<std::string> fileNames;
        std::vector<std::string> lc_betaIntNames;
        std::vector<std::string> betaIntNames;

        // Indices
        std::unordered_map<std::string, int> nheader;
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
        std::unordered_map<std::string, std::vector<int>> betaIntColumn2;
        std::unordered_map<std::string, std::vector<int>> mb_covIntColumn;
        std::unordered_map<std::string, std::vector<int>> rb_covIntColumn;      
};

void processFileHeader(int nInt1, int nInt2, bool mb, bool rb, std::vector<std::string> lc_intNames, std::vector<std::string> lc_intNames2, std::vector<std::string> fileNames, FileInfo* fip);
void printOutputHeader(bool mb, bool rb, std::string output, size_t nInt1, std::vector<std::string> intNames); 
void printHeaderMissingError(std::string fileName, std::string column);

#endif
