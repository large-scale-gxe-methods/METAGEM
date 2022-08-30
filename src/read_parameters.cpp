/****************************************************************************
  ReadParameters.cpp reads parameters from a file
****************************************************************************/
#include "metagem.h"
#include <boost/program_options.hpp>

#define VERSION "1.0"

namespace po = boost::program_options;

void print_help();

// Function to process command line arguments
void CommandLine::processCommandLine(int argc, char* argv[]) {


    cout << "\n*********************************************************\n";
    cout << "Welcome to METAGEM v" << VERSION << "\n";
    cout << "(C) 2021 Duy Pham and Han Chen \n";
    cout << "GNU General Public License v3\n";
    cout << "*********************************************************\n";


    // GEM parameters. Details are printed from the print_help() function below.

    // General commands
    po::options_description general("General options");
    general.add_options()
        ("help, h", "")
        ("version", "");

    // Input/Output file options
    po::options_description files("Input/Output file options");
    files.add_options()
        ("meta-files", po::value<std::vector<std::string>>()->multitoken(), "")
        ("exposure-names", po::value<std::vector<std::string>>()->multitoken(), "")
        ("int-covar-names", po::value<std::vector<std::string>>()->multitoken(), "")
        ("out", po::value<std::string>()->default_value("metagem.out"), "");

    // Combine all options together
    po::options_description all("Options");
    all.add(general).add(files);

    po::variables_map out;


    // QC
    try {
        po::store(po::command_line_parser(argc, argv)
            .options(all)
            .style(po::command_line_style::unix_style
                | po::command_line_style::allow_long_disguise)
            .run(),
            out);

    }
    catch (po::error const& e) {
        std::cerr << e.what() << endl;
        exit(EXIT_FAILURE);
    }
    po::notify(out);



    // General
    if (out.count("help")) {
        print_help();
        exit(1);

    }
    if (out.count("version")) {
        cout << "\nMETAGEM version: " << VERSION << "\n\n.";
        exit(1);

    }


    // Input/Output Files
    if (out.count("meta-files")) {
        fileNames = out["meta-files"].as<std::vector<std::string>>();

        size_t nfiles = fileNames.size();
        if (nfiles < 2) {
            cerr << "\nERROR: METAGEM requires at least 2 GEM results file.\n\n";
            exit(1);
        }

        std::set<std::string> s(fileNames.begin(), fileNames.end());
        if (s.size() != nfiles) 
        {
            cerr << "\nERROR: There are duplicate file names (--files).\n\n";
            exit(1);
        }

    }
    else {
        cerr << "\nERROR: Meta files (--meta-files) are needed. \n\n";
        exit(1);

    }

    // Output file
    if (out.count("out")) {
        outFile = out["out"].as<std::string>();

        std::ofstream results(outFile, std::ofstream::binary);
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
            boost::filesystem::remove(outFile.c_str());
            exit(1);
        }
        results.close();

        boost::filesystem::remove(outFile.c_str());
    }
    

    // Exposures
    if (out.count("exposure-names")) {
        intNames = out["exposure-names"].as<std::vector<std::string>>();
        nExp = intNames.size();

        std::set<std::string> s(intNames.begin(), intNames.end());
        if (s.size() != intNames.size()) 
        {
            cerr << "\nERROR: There are duplicate exposure names (--exposure-names).\n\n";
            exit(1);
        }
    }

    // Interaction covariates
    if (out.count("int-covar-names")) {
        std::vector<std::string> icovNames = out["int-covar-names"].as<std::vector<std::string>>();
        nCov = icovNames.size();

        std::set<std::string> s(icovNames.begin(), icovNames.end());
        if (s.size() != nCov) 
        {
            cerr << "\nERROR: There are duplicate interation covariate names (--int-covar-names).\n\n";
            exit(1);
        }

        intNames.insert(intNames.end(), icovNames.begin(), icovNames.end());
    }

}






void print_help() {

 
    cout << "General Options: " << endl
        << "   --help \t\t Prints available options and exits." << endl
        << "   --version \t\t Prints the version of METAGEM and exits." << endl;
    cout << endl << endl;



    cout << "Input File Options: " << endl
        << "   --meta-files \t\t Output files from GEM 'meta' or 'full' option for meta-analysis." << endl
        << "   --exposure-names \t " << endl
        << "   --int-covar-names \t" << endl
        << "   --out \t\t Full path and extension to where METAGEM output results. \n \t\t\t    Default: metagem.out" << endl;
    cout << endl << endl;




    cout << endl << endl;

}