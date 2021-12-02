#include "metagem.h"

void processFileHeader(int Sq1, std::vector<std::string> lc_intNames, std::string fileName, FileInfo* fip) 
{

	// Read input file
	std::ifstream file;
	file.open(fileName);
	if (!file.is_open()) {
		cerr << "\nERROR: Cannot open the file: " << fileName << "\n\n";
		exit(1);
	}
	fip->fileNames.push_back(fileName);

	// Ignore comments in the text file (ex. #dispersion: )
	std::string line;
	getline(file, line);
	while(line.rfind("#", 0) == 0)
	{
		getline(file, line);
	}

	// Check if there are any results in the file
	int nvars = 0;
	while (getline(file, line)) nvars++;

	if (nvars == 0 ) {
		cerr << "\nERROR: No results in the file: " << fileName << ".\n\n";
		file.close();
		exit(1);
	} else {
		file.clear(); 
		file.seekg(0, file.beg);
		getline(file, line);
		while(line.rfind("#", 0) == 0)
		{
			getline(file, line);
		}
	}
	fip->nvars[fileName] = nvars;

	
	// Read the column headers
	std::unordered_map<std::string, int> columnNames;

	int header_i = 0;
	std::string header;
	std::istringstream iss(line);
	while (getline(iss, header, '\t')) 
	{
		header.erase(std::remove(header.begin(), header.end(), '\r'), header.end());

		if (columnNames.find(header) != columnNames.end()) {
			cerr << "\nERROR: There are duplicate header names (" << header << ") in file [" << fileName << "].\n\n";
			file.close();
			exit(1);
		}
		columnNames[header] = header_i;

		header_i++;
	}
	file.close();

	
	// Get variant information columns
	if (columnNames.find("SNPID") != columnNames.end()) 
	{
		fip->snpColumn[fileName] = columnNames["SNPID"];
	} 
	else 
	{
		cerr << "\nERROR: The file [" << fileName << "] does not contain a SNPID column.\n\n";
		exit(1);
	}

	if (columnNames.find("CHR") != columnNames.end())
	{
		fip->chrColumn[fileName] = columnNames["CHR"];
	}
	else
	{
		cerr << "\nERROR: The file [" << fileName << "] does not contain a CHR column.\n\n";
		exit(1);
	}

	if (columnNames.find("POS") != columnNames.end())
	{
		fip->posColumn[fileName] = columnNames["POS"];
	}
	else
	{
		cerr << "\nERROR: The file [" << fileName << "] does not contain a POS column.\n\n";
		exit(1);
	}

	if (columnNames.find("Non_Effect_Allele") != columnNames.end())
	{
		fip->nonEffectColumn[fileName] = columnNames["Non_Effect_Allele"];
	}
	else
	{
		cerr << "\nERROR: The file [" << fileName << "] does not contain a Non_Effect_Allele column.\n\n";
		exit(1);
	}

	if (columnNames.find("Effect_Allele") != columnNames.end())
	{
		fip->effectColumn[fileName] = columnNames["Effect_Allele"];
	}
	else
	{
		cerr << "\nERROR: The file [" << fileName << "] does not contain a Effect_Allele column.\n\n";
		exit(1);
	}

	if (columnNames.find("N_Samples") != columnNames.end())
	{
		fip->nSampleColumn[fileName] = columnNames["N_Samples"];
	}
	else
	{
		cerr << "\nERROR: The file [" << fileName << "] does not contain a N_Samples column.\n\n";
		exit(1);
	}

	if (columnNames.find("AF") != columnNames.end())
	{
		fip->freqColumn[fileName] = columnNames["AF"];
	}
	else
	{
		cerr << "\nERROR: The file [" << fileName << "] does not contain a AF column.\n\n";
		exit(1);
	}
	

	// Get the Beta Marginal column
	if (columnNames.find("Beta_Marginal") != columnNames.end())
	{
		fip->betaMargColumn[fileName] = columnNames["Beta_Marginal"];
	}
	else
	{
		cerr << "\nERROR: The file [" << fileName << "] does not contain a Beta_Marginal column.\n\n";
		exit(1);
	}

	// Get marginal model based SE columns if it exists
	if (columnNames.find("SE_Beta_Marginal") != columnNames.end())
	{
		fip->mb_seMargColumn[fileName] = columnNames["SE_Beta_Marginal"];
	}
	else 
	{
		cerr << "\nERROR: The file [" << fileName << "] does not contain model-based marginal SE.\n\n";
		exit(1);
	}

	// Get marginal robust SE columns if it exists
	if (columnNames.find("robust_SE_Beta_Marginal") != columnNames.end())
	{
		fip->rb_seMargColumn[fileName] = columnNames["robust_SE_Beta_Marginal"];
	}
	else 
	{
		cerr << "\nERROR: The file [" << fileName << "] does not contain robust marginal SE.\n\n";
		exit(1);
	}

		

	// Get Beta Interaction columns
	std::vector<std::string> betaIntNames;

	if (columnNames.find("Beta_G") == columnNames.end()) {
		cerr << "\nERROR: The file [" << fileName << "] does not contain a Beta_G column.\n\n";
		exit(1);
	}
	betaIntNames.push_back("G");

	int nints = 1;
	for (std::pair<std::string, int> element : columnNames)
	{
		std::string key = element.first;
		if (key.rfind("Beta_G-", 0) == 0) {
			key.erase(0, 5);
			betaIntNames.push_back(key);
			nints++;
		}
	}
	if (nints <= 1 ) {
		cerr << "\nERROR: The file [" << fileName << "] does not contain any Beta_G-* columns.\n\n";
		exit(1);
	}

	// Convert to lowercase strings for matching
	std::vector<std::string> lc_betaIntNames = betaIntNames;
	for(std::string &s : lc_betaIntNames){
		std::transform(s.begin(), s.end(), s.begin(), [](char c){ return std::tolower(c); });
	}
	

	if(nints != Sq1)
	{
		cerr << "\nERROR: The file [" << fileName << "] does not contain the same number of interaction terms.\n\n";
		exit(1);
	}
		
	// Order the interactions base on the first set of interaction vector
	std::vector<int> betaIntColumn;
	std::vector<std::string> ord_betaIntNames;
	for (int i = 0; i < Sq1; i++)
	{
		auto it = std::find(lc_betaIntNames.begin(), lc_betaIntNames.end(), lc_intNames[i]);
		if (it != lc_betaIntNames.end())
		{
			auto idx = std::distance(lc_betaIntNames.begin(), it);
			ord_betaIntNames.push_back(betaIntNames[idx]);
			betaIntColumn.push_back(columnNames["Beta_" + betaIntNames[idx]]);
		}
		else 
		{
			cerr << "\nERROR: The file [" << fileName << "] does not contain the GxE term: " << fip->betaIntNames[i] << ".\n\n";
			exit(1);
		}
	}
	fip->betaIntColumn[fileName] = betaIntColumn;

	
	// Get columns containing the model-based summary statistics for the interaction terms	
	std::vector<int> mb_covIntColumn(nints * nints);
	for (int i = 0; i < nints; i++) 
	{
		std::string s1 = "SE_Beta_" + ord_betaIntNames[i];
		if (columnNames.find(s1) != columnNames.end()) 
		{
			mb_covIntColumn[i*nints + i] = columnNames[s1];
		}
		else 
		{
			cerr << "\nERROR: The file [" << fileName << "] does not contain the column " << s1 << ".\n\n";
			exit(1);
		}
	}

	for (int i = 0; i < nints; i++) 
	{
		for (int j = i+1; j < nints; j++) 
		{
			int idx;
			std::string s2 = "Cov_Beta_" + ord_betaIntNames[i] + "_" + ord_betaIntNames[j];
			std::string s3 = "Cov_Beta_" + ord_betaIntNames[j] + "_" + ord_betaIntNames[i];
			if (columnNames.find(s2) != columnNames.end()) 
			{
				idx = columnNames[s2];
			} 
			else if (columnNames.find(s3) != columnNames.end())
			{
				idx = columnNames[s3];
			}
			else 
			{
				cerr << "\nERROR: The file [" <<  fileName << "] does not contain the column " << s2 << " or " << s3 << ".\n\n";
				exit(1);
			}
			mb_covIntColumn[(i * nints) + j] = idx;
			mb_covIntColumn[(j * nints) + i] = idx;
		}
	}
	fip->mb_covIntColumn[fileName] = mb_covIntColumn;

	
	// Get columns containing the robust summary statistics for the interaction terms	
	std::vector<int> rb_covIntColumn(nints * nints);
	for (int i = 0; i < nints; i++) 
	{
		std::string s1 = "robust_SE_Beta_" + ord_betaIntNames[i];
		if (columnNames.find(s1) != columnNames.end()) 
		{
			rb_covIntColumn[i*nints + i] = columnNames[s1];
		}
		else 
		{
			cerr << "\nERROR: The file [" << fileName << "] does not contain the column " << s1 << ".\n\n";
			exit(1);
		}
	}

	for (int i = 0; i < nints; i++) 
	{
		for (int j = i+1; j < nints; j++) 
		{
			int idx;
			std::string s2 = "robust_Cov_Beta_" + ord_betaIntNames[i] + "_" + ord_betaIntNames[j];
			std::string s3 = "robust_Cov_Beta_" + ord_betaIntNames[j] + "_" + ord_betaIntNames[i];

			if (columnNames.find(s2) != columnNames.end()) 
			{
				idx = columnNames[s2];
			} 
			else if (columnNames.find(s3) != columnNames.end())
			{
				idx = columnNames[s3];
			}
			else 
			{
				cerr << "\nERROR: The file [" <<  fileName << "] does not contain the column " << s2 << " or " << s3 << ".\n\n";
				exit(1);
			}
			rb_covIntColumn[(i * nints) + j] = idx;
			rb_covIntColumn[(j * nints) + i] = idx;
		}
	}
	fip->rb_covIntColumn[fileName] = rb_covIntColumn;


	cout << "Read file [" << fileName << "] with " << std::to_string(nvars) << " variants.\n";
}



void printOutputHeader(std::string output, size_t Sq1, std::vector<std::string> intNames) 
{
	for(std::string &s : intNames){
		s = "G-" + s;
	}
	intNames.insert(intNames.begin(), "G");


	std::ofstream results(output, std::ofstream::binary);

	results << "SNPID" << "\t" << "CHR" << "\t" << "POS" << "\t" << "Non_Effect_Allele" << "\t" << "Effect_Allele" << "\t" << "N_Samples" << "\t" << "AF" << "\t";
	results << "Beta_Marginal" << "\t" << "SE_Beta_Marginal" << "\t" << "robust_SE_Beta_Marginal" << "\t";

	// Print beta header
	for (size_t i = 0; i < Sq1; i++) 
	{
		results << "Beta_" << intNames[i] << "\t";
	}

	// Print model-based covariance
	for (size_t i = 0; i < Sq1; i++) 
	{
		results << "SE_Beta_" << intNames[i] << "\t"; 
	}

	for (size_t i = 0; i < Sq1; i++) 
	{
		for (size_t j = 0; j < Sq1; j++) 
		{
			if (i < j) 
			{
				results << "Cov_Beta_" << intNames[i] << "_" << intNames[j] << "\t"; 
			}
		}
	}

	// Print robust covariance
	for (size_t i = 0; i < Sq1; i++) 
	{
		results << "robust_SE_Beta_" << intNames[i] << "\t"; 
	}

	for (size_t i = 0; i < Sq1; i++) 
	{
		for (size_t j = 0; j < Sq1; j++) 
		{
			if (i < j) 
			{
				results << "robust_Cov_Beta_" << intNames[i] << "_" << intNames[j] << "\t"; 
			}
		}
	}

	// Print p-value
	results << "P_value_Marginal" << "\t" << "P_Value_Interaction" << "\t" << "P_Value_Joint" << "\t";
	results << "robust_P_value_Marginal" << "\t" << "robust_P_Value_Interaction" << "\t" << "robust_P_Value_Joint\n";
	
	results.close();
}
