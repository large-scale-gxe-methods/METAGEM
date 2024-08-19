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
        std::unordered_map<std::string, std::vector<int>> betaIntColumn3;
        std::unordered_map<std::string, std::vector<int>> mb_covIntColumn;
        std::unordered_map<std::string, std::vector<int>> mb_covIntColumn2;
        std::unordered_map<std::string, std::vector<int>> mb_covIntColumn3;
        std::unordered_map<std::string, std::vector<int>> rb_covIntColumn;      
        std::unordered_map<std::string, std::vector<int>> rb_covIntColumn2;
        std::unordered_map<std::string, std::vector<int>> rb_covIntColumn3;
};

void processFileHeader(int nInt1, int nInt2, int nint3, bool mb, bool rb, bool additionalJoint, bool additionalInteraction, bool renameHeaders, std::vector<std::string> lc_intNames, std::vector<std::string> lc_intNames2, std::vector<std::string> lc_intNames3, std::vector<std::string> fileNames, std::map<std::string, std::map<std::string, std::string>> headerRenamings, FileInfo* fip);
void printOutputHeader(bool mb, bool rb, bool additionalJoint, bool additionalInteraction, std::string output, std::string output2, std::string output3, size_t nInt1, size_t nInt2, size_t nInt3, std::vector<std::string> intNames, std::vector<std::string> intNames2, std::vector<std::string> intNames3); 
void printHeaderMissingError(std::string fileName, std::string column);

#endif
