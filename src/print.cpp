#include "metagem.h"


void printWelcome() {
    cout << "\n*********************************************************\n";
    cout << "Welcome to METAGEM v" << VERSION << "\n";
    cout << "(C) 2021-2022 Duy Pham and Han Chen \n";
    cout << "GNU General Public License v3\n";
    cout << "*********************************************************\n";
}

void printHeaderMissingError(std::string fileName, std::string column) {
    cerr << "\nERROR: The file [" << fileName << "] does not contain a " << column << " column.\n\n";
}

void printNColumnError(int z, int n, std::string fileName, int nheader) {
  cerr << "\nERROR: Unexpect number of columns (" << z << ") for variant number " << n+1 << "in file [" << fileName << "]. Expected: " << nheader << ".\n\n";
}

void printZeroVariantsError(std::string fileName) {
  cerr << "\nERROR: No variants in file [" << fileName << "].\n\n";
}

void printOpenFileError(std::string fileName) {
    cerr << "\nERROR: Cannot open the file [" << fileName << "].\n\n";
}

void printProcessingFiles() {
  cout << "Processing files...\n";
}

void printProcessedFiles(std::string fileName, int n) {
  cout << "Processed file [" << fileName << "] with " << n << " variants.\n";
}

void printMetaBegin(int nFiles, int nvars) {
  cout << "Performing meta-analysis with " << std::to_string(nFiles) << " files and " << std::to_string(nvars) << " variants...\n";
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

void printTimeCompleted(double wall0, double wall1, double cpu0, double cpu1) {
    cout << "\n*********************************************************\n";
    cout << "Wall Time = " << wall1 - wall0 << " (sec)\n";
    cout << "CPU Time  = " << cpu1  - cpu0  << " (sec)\n\n";
}