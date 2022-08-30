/****************************************************************************
  ReadParameters.cpp reads parameters from a file
****************************************************************************/
#include "metagem.h"


void print_help();

// Function to process command line arguments
void CommandLine::processCommandLine(int argc, char* argv[]) {


    cout << "\n*********************************************************\n";
    cout << "Welcome to METAGEM v" << VERSION << "\n";
    cout << "(C) 2021 Duy Pham and Han Chen \n";
    cout << "GNU General Public License v3\n";
    cout << "*********************************************************\n";


    // GEM options. Details are printed from the print_help() function below.
    CLI::App app{""};

    // Defaults
    int metaOpt_in = 0;
    std::string outFile_in = "metagem.out";

    app.add_option("--input-files", fileNames, "")->expected(0, 1000000);
    app.add_option("--input-file-list", metaFileList, "")->expected(1);
    app.add_option("--exposure-names", intNames, "")->expected(1, 1000000)->required();
    app.add_option("--out", outFile_in, "")->expected(1);
    app.add_option("--meta-option", metaOpt_in, "");

    try
    {
        app.parse( argc, argv);

        // Input files
        size_t fns = fileNames.size();
        size_t fls = metaFileList.length();
        if (fns == 0 && fls == 0) {
            cerr << "\nERROR: --input-files or --input-file-list is required.\n\n";
            exit(1);

        } else if (fns > 0 && fls > 0) {
            cerr << "\nERROR: Both --input-files and --input-file-list are specified.\n\n";
            exit(1);

        } else if (fns > 0) {
            if (fns < 2) {
                cerr << "\nERROR: METAGEM requires at least 2 input files.\n\n";
                exit(1);                
            }

            std::set<std::string> s(fileNames.begin(), fileNames.end());
            if (s.size() != fns) 
            {
                cerr << "\nERROR: There are duplicate input file names.\n\n";
                exit(1);
            }

        } else if (fls > 0) {
            std::ifstream file;
            std::string line;
            file.open(metaFileList);
            if (!file.is_open()) {
                cerr << "\nERROR: Cannot open the file: [" << metaFileList << "].\n\n";
                exit(1);
            }

            size_t nfiles = 0;
            while(getline(file, line)) {
                line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
                fileNames.push_back(line);
                nfiles++;
            }

            if (nfiles <= 2)
            {
                cerr << "\nERROR: METAGEM requires at least 2 input files.\n\n";
                exit(1);
            }
        
            std::set<std::string> s(fileNames.begin(), fileNames.end());
            if (s.size() != nfiles) 
            {
                cerr << "\nERROR: There are duplicate input file names.\n\n";
                exit(1);
            }

        } else {
            cerr << "\nERROR: Unrecognized combination of --input-files and --input-file-list.\n\n";
            exit(1);
        }


        // Exposure names
        std::set<std::string> s(intNames.begin(), intNames.end());
        if (s.size() != intNames.size()) {
            cerr << "\nERROR: There are duplicate exposure names (--exposure-names).\n\n";
            exit(1);
        }
        nExp = intNames.size();


        // Output file
        outFile = outFile_in;
        std::ofstream results(outFile);
        if (!results) {
            cerr << "\nERROR: Output file could not be opened.\n\n";
            exit(1);
        }

        if (results.fail()) {
            cerr << "\nERROR: Output file could not be opened.\n\n";
            exit(1);
        }

        results << "test" << endl;
        if (results.fail()) {
            cerr << "\nERROR: Cannot write to output file.\n\n";
            results.close();
            
            if (std::remove(outFile.c_str()) != 0) {
                cerr << "\nERROR: Cannot delete output file.\n\n";
            }
            exit(1);
        }
        results.close();
        
        if (std::remove(outFile.c_str()) != 0) {
            cerr << "\nERROR: Cannot delete output file.\n\n";
            exit(1);
        }

        
        // Meta Option
        metaOpt = metaOpt_in;
        if (metaOpt < 0 || metaOpt > 2) {
            cerr << "\nERROR: The --meta-option integer value must be 0, 1, or 2.\n\n";
            exit(1);
        }

        if (metaOpt == 0) {
            mb = true;
            rb = true;

        } else if (metaOpt == 1) {
            mb = true;
        } else if (metaOpt == 2) {
            rb = true;
        } else {
            cerr << "\nERROR: The --meta-option integer value must be 0, 1, or 2.\n\n";
            exit(1);
        }
        
    }
    catch( const CLI::CallForHelp &e )
    {
        print_help();
        exit(1);
    }
}





void print_help() {

 
    cout << "General Options: " << endl
        << "   --help \t\t Prints available options and exits." << endl;
    cout << endl << endl;



    cout << "Input File Options: " << endl
        << "   --input-files \t Output files from GEM 'meta' or 'full' option." << endl
        << "   --input-file-list \t A no header text file containing a single file name per line." << endl
        << "   --exposure-names \t The names of the exposures." << endl
        << "   --out \t\t Full path and extension to where METAGEM output results. \n \t\t\t    Default: metagem.out" << endl
        << "   --meta-option \t Integer value indicating which summary statistics should be used for meta-analysis. \n\t\t\t    0: Both model-based and robust summary statistics. \n \t\t\t    1: model-based summary statistics. \n \t\t\t    2: robust summary statistics. \n \t\t\t    Default: 0" << endl;
    cout << endl << endl;




    cout << endl << endl;

}