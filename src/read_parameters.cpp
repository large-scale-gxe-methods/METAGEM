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
        ("input-files", po::value<std::vector<std::string>>()->multitoken(), "")
        ("input-file-list", po::value<std::string>(), "")       
        ("exposure-names", po::value<std::vector<std::string>>()->multitoken(), "")
        ("int-covar-names", po::value<std::vector<std::string>>()->multitoken(), "")
        ("meta-option", po::value<int>()->default_value(0))
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
    if (out.count("input-files")) {
        fileNames = out["input-files"].as<std::vector<std::string>>();

        size_t nfiles = fileNames.size();
        if (nfiles < 2) {
            cerr << "\nERROR: METAGEM requires at least 2 input files.\n\n";
            exit(1);
        }

        std::set<std::string> s(fileNames.begin(), fileNames.end());
        if (s.size() != nfiles) 
        {
            cerr << "\nERROR: There are duplicate input file names.\n\n";
            exit(1);
        }
                for (int s = 0; s < nfiles; s++)
        {
            cout << fileNames[s] << endl;
        }

    }
    if (out.count("input-file-list")) {
        if (out.count("input-files")) {
            cerr << "\nERROR: --input-files is also specified.\n\n";
            exit(1);
        }
        metaFileList = out["input-file-list"].as<std::string>();

        // Read input file
		std::ifstream file;
        std::string line;
		file.open(metaFileList);
		if (!file.is_open()) {
			cerr << "\nERROR: Cannot open the file: " << metaFileList << "\n\n";
			exit(1);
		}

        size_t nfiles = 0;
		while(getline(file, line))
		{
			line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
            fileNames.push_back(line);
            nfiles++;
		}

        if (nfiles <= 0)
        {
            cerr << "\nERROR: No file names in --input-file-list.\n\n";
            exit(1);
        }
        
        if (nfiles < 2)
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
        for (int s = 0; s < nfiles; s++)
        {
            cout << fileNames[s] << endl;
        }
    }
    if (!out.count("input-files") && !out.count("input-file-list")) {
        cerr << "\nERROR: --input-files or --input-file-list is required.\n\n";
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


    if (out.count("meta-option")) {
        metaOpt = out["meta-option"].as<int>();

        if (metaOpt < 0 || metaOpt > 2)
        {
            cerr << "\nERROR: The --meta-option should be 0, 1, or 2.\n\n";
            exit(1);
        }

        if (metaOpt == 0)
        {
            mb = true;
            rb = true;
        }
        else if (metaOpt == 1)
        {
            mb = true;
        }
        else if (metaOpt == 2)
        {
            rb = true;
        }
    }
}






void print_help() {

 
    cout << "General Options: " << endl
        << "   --help \t\t Prints available options and exits." << endl
        << "   --version \t\t Prints the version of METAGEM and exits." << endl;
    cout << endl << endl;



    cout << "Input File Options: " << endl
        << "   --input-files \t Output files from GEM 'meta' or 'full' option." << endl
        << "   --input-file-list \t A no header text file containing a single file name per line." << endl
        << "   --exposure-names \t " << endl
        << "   --out \t\t Full path and extension to where METAGEM output results. \n \t\t\t    Default: metagem.out" << endl
        << "   --meta-option \t Integer value indicating which summary statistics should be used for meta-analysis. \n\t\t\t    0: Both model-based and robust summary statistics. \n \t\t\t    1: model-based summary statistics. \n \t\t\t    2: robust summary statistics. \n \t\t\t    Default: 0" << endl;
    cout << endl << endl;




    cout << endl << endl;

}