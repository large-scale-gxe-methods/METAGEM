#include "metagem.h"

int main(int argc, char* argv[]) 
{
	// Process command line arguments
	CommandLine cmd;
	cmd.processCommandLine(argc, argv);

	metagem(cmd);

	return(0);
}


int flipAllele(std::string a1, std::string a2, std::string b1, std::string b2)
{
	std::transform(b1.begin(), b1.end(), b1.begin(), [](char c){ return std::tolower(c); });
	std::transform(b2.begin(), b2.end(), b2.begin(), [](char c){ return std::tolower(c); });

	if ((a1 == b1) && (a2 == b2))
	{
		return 0;
	}
	if ((a1 == b2) && (a2 == b1))
	{
		return 1;
	}

	return 2;
}



void metagem(CommandLine cmd)
{
	size_t nFiles = cmd.fileNames.size();
	size_t nExp   = cmd.nExp;
	size_t nExp1  = nExp + 1; // plus G


	int Sq1 = 1;
	std::vector<std::string> lc_intNames = cmd.intNames;
	for(std::string &s : lc_intNames){
		std::transform(s.begin(), s.end(), s.begin(), [](char c){ return std::tolower(c); });
		s = "g-" + s;
		Sq1++;
	}
	lc_intNames.insert(lc_intNames.begin(), "g");



	// Read the header of each file
	FileInfo* fip = new FileInfo();
	for (size_t i = 0; i < nFiles; i++) 
	{
		processFileHeader(Sq1, lc_intNames, cmd.fileNames[i], fip);
	}
	
	cout << "\nPerforming meta-analysis with " << std::to_string(nFiles) << " files.\n";

	int nvars = fip->nvars[fip->fileNames[0]];
	int uSize = nvars * Sq1;
	int vSize = nvars * (Sq1 * Sq1);

	std::unordered_map<std::string, int> snpid_idx;
	std::vector<std::string> snpid(nvars);
	std::vector<std::string> chr(nvars);
	std::vector<std::string> pos(nvars);
	std::vector<std::string> effectAllele(nvars);
	std::vector<std::string> nonEffectAllele(nvars);
	std::vector<std::string> lc_effectAllele(nvars);
	std::vector<std::string> lc_nonEffectAllele(nvars);
	std::vector<double> nSamples(nvars);
	std::vector<double> AF(nvars);
	
	std::vector<int> skip(nvars);
	std::vector<size_t> fileCount(nvars);

	// Get U vectors and V matrices from each file
	std::vector<double> fileBetaMarg(nvars);
	std::vector<double> file_mbMV(nvars);
	std::vector<double> file_rbMV(nvars);

	std::vector<double> fileBetaInt(uSize);
	std::vector<double> file_mbV(vSize);
	std::vector<double> file_rbV(vSize);
	double* file_mbU  = new double[uSize];
	double* file_rbU  = new double[uSize];

	double* mb_MU = new double[uSize];
	double* rb_MU = new double[uSize];
	double* mb_MV = new double[uSize];
	double* rb_MV = new double[uSize];
	double* mb_U  = new double[uSize];
	double* rb_U  = new double[uSize];
	double* mb_V  = new double[vSize];
	double* rb_V  = new double[vSize];

	boost::math::chi_squared chisq_dist_M(1);
	boost::math::chi_squared chisq_dist_Int(nExp);
	boost::math::chi_squared chisq_dist_Joint(nExp1);

	for (size_t f = 0; f < nFiles; f++)
	{

		std::string fileName = fip->fileNames[f];
	
		int nrow            = fip->nvars[fileName];
		int snpColumn       = fip->snpColumn[fileName];
		int chrColumn       = fip->chrColumn[fileName];
		int posColumn       = fip->posColumn[fileName];
		int effectColumn    = fip->effectColumn[fileName];
		int nonEffectColumn = fip->nonEffectColumn[fileName];
		int nSampleColumn   = fip->nSampleColumn[fileName];
		int afColumn        = fip->freqColumn[fileName];
		int betaMargColumn  = fip->betaMargColumn[fileName];
		int mb_seMargColumn = fip->mb_seMargColumn[fileName];
		int rb_seMargColumn = fip->rb_seMargColumn[fileName];		

		std::vector<int> betaIntColumn   = fip->betaIntColumn[fileName];
		std::vector<int> mb_covIntColumn = fip->mb_covIntColumn[fileName];
		std::vector<int> rb_covIntColumn = fip->rb_covIntColumn[fileName];

		// Read input file
		std::ifstream file;
		file.open(fileName);
		if (!file.is_open()) {
			cerr << "\nERROR: Cannot open the file: " << fileName << "\n\n";
			exit(1);
		}

		// Ignore comments in the text file (ex. #dispersion: )
		std::string line;
		getline(file, line);
		while(line.rfind("#", 0) == 0)
		{
			getline(file, line);
		}


		int n  = 0;
		int bi = 0;
		int ci = 0;
		if (f == 0)
		{
			while(getline(file, line))
			{
				std::istringstream iss(line);

				std::string value;
				std::vector <std::string> values;
				while (getline(iss, value, '\t')) values.push_back(value);


				snpid[n] = values[snpColumn];
				chr[n]   = values[chrColumn];
				pos[n]   = values[posColumn];

				nSamples[n] = std::stod(values[nSampleColumn]);
 				AF[n]       = std::stod(values[afColumn]) * 2.0 * nSamples[n];

				std::string a1 = values[nonEffectColumn];
				std::string a2 = values[effectColumn];
				nonEffectAllele[n] = a1;
				effectAllele[n]    = a2;

				std::transform(a1.begin(), a1.end(), a1.begin(), [](char c){ return std::tolower(c); });
				std::transform(a2.begin(), a2.end(), a2.begin(), [](char c){ return std::tolower(c); });
				lc_nonEffectAllele[n] = a1;
				lc_effectAllele[n]    = a2;

				fileBetaMarg[n] = std::stod(values[betaMargColumn]);
				file_mbMV[n]    = std::stod(values[mb_seMargColumn]);
				file_rbMV[n]    = std::stod(values[rb_seMargColumn]);

				for (int i = 0; i < Sq1; i++) {
					fileBetaInt[bi] = std::stod(values[betaIntColumn[i]]);
					bi++;
				}

				int idx = n * (Sq1 * Sq1);
				for (int i = 0; i < Sq1; i++) {
					int ii = i * Sq1;
					int ni = idx + ii + i;
					for (int j = 0; j < Sq1; j++) {
						file_mbV[ci] = std::stod(values[mb_covIntColumn[ii + j]]);
						file_rbV[ci] = std::stod(values[rb_covIntColumn[ii + j]]);
						ci++;
					}
					file_mbV[ni] *= file_mbV[ni];
					file_rbV[ni] *= file_rbV[ni];
				}

				snpid_idx[snpid[n]] = n;
				n++;
			}
		}
		else 
		{
			while(getline(file, line))
			{
				std::string value;
				std::vector <std::string> values;
				std::istringstream iss(line);
				while (getline(iss, value, '\t')) values.push_back(value);

				int val;
				double dir = 1.0;
				std::unordered_map<std::string, int>::iterator it = snpid_idx.find(values[snpColumn]);
				if (it != snpid_idx.end())
				{
					val = it->second;
					int flip = flipAllele(lc_nonEffectAllele[val], lc_effectAllele[val], values[nonEffectColumn], values[effectColumn]);
					double ns = std::stod(values[nSampleColumn]);
					nSamples[val] += ns;				
					if (flip == 0)
					{
						AF[val] += std::stod(values[afColumn]) * 2.0 * ns;
					}
					else if (flip == 1)
					{
						dir = -1.0;
						AF[val] += (1 - std::stod(values[afColumn])) * 2.0 * ns;
					}
					else
					{
						continue;
					}
					fileCount[val]++;
				}
				else {
					continue;
				}

				fileBetaMarg[val] = std::stod(values[betaMargColumn]);
				file_mbMV[val]    = std::stod(values[mb_seMargColumn]);
				file_rbMV[val]    = std::stod(values[rb_seMargColumn]);

				for (int i = 0; i < Sq1; i++) {
					fileBetaInt[bi] = std::stod(values[betaIntColumn[i]]) * dir;
					bi++;
				}

				int idx = val * (Sq1 * Sq1);
				for (int i = 0; i < Sq1; i++) {
					int ii = i * Sq1;
					int ni = idx + ii + i;
					for (int j = 0; j < Sq1; j++) {
						file_mbV[ci] = std::stod(values[mb_covIntColumn[ii + j]]);
						file_rbV[ci] = std::stod(values[rb_covIntColumn[ii + j]]);
						ci++;
					}
					file_mbV[ni] *= file_mbV[ni];
					file_rbV[ni] *= file_rbV[ni];
				}
				n++;
			}

		}
		file.close();

		if (n != nrow) {
			cerr << "\nERROR: The number of variants read in [" << fileName << "] does not match original count.\n\n";
			exit(1);
		}

		for (int i = 0; i < nvars; i++)
		{
			if (fileCount[i] != f)
			{
				skip[i] = 1;
				continue;
			}

			int is  = (i * Sq1);
			int iss = (i * Sq1 * Sq1);

			subMatInv(&file_mbV[0], Sq1, iss);
			subMatInv(&file_rbV[0], Sq1, iss);
			subMatsubVecprod(&file_mbV[0], &fileBetaInt[0], file_mbU, Sq1, Sq1, iss, is, is);
			subMatsubVecprod(&file_rbV[0], &fileBetaInt[0], file_rbU, Sq1, Sq1, iss, is, is);

			// Marginal
			file_mbMV[i] *= file_mbMV[i];
			file_rbMV[i] *= file_rbMV[i];
			file_mbMV[i] = 1 / file_mbMV[i];
			file_rbMV[i] = 1 / file_rbMV[i];

			mb_MU[i] += file_mbMV[i] * fileBetaMarg[i];
			rb_MU[i] += file_rbMV[i] * fileBetaMarg[i];
			mb_MV[i] += file_mbMV[i];
			rb_MV[i] += file_rbMV[i];
		}

		matAdd(mb_U, file_mbU, uSize, 1);
		matAdd(rb_U, file_rbU, uSize, 1);
		matAdd(mb_V, &file_mbV[0], vSize, 1);
		matAdd(rb_V, &file_rbV[0], vSize, 1);
	}
	file_mbV.clear();
	file_rbV.clear();
	file_mbMV.clear();
	file_rbMV.clear();
	fileBetaInt.clear();
	delete[] file_mbU;
	delete[] file_rbU;


	// Create the output file and write column header names
	printOutputHeader(cmd.outFile, Sq1, cmd.intNames);
	std::ofstream results(cmd.outFile, std::ios_base::app);
	std::ostringstream oss;

	double* betaInt = new double[Sq1];
	double* StempE  = new double[nExp];
	double* StempGE = new double[nExp1];
	double* Ai      = new double[nExp1 * nExp1];
	double* VE      = new double[nExp * nExp];
	for (int i = 0; i < nvars; i++)
	{
		if (skip[i])
		{
			continue;
		}

		int is  = (i * Sq1);
		int iss = (i * Sq1 * Sq1);

		AF[i] = AF[i] / 2.0 / nSamples[i];

		subMatInv(mb_V, Sq1, iss);
		subMatInv(rb_V, Sq1, iss);

		// Marginal Effect
		double betaMarg = mb_MU[i] / mb_MV[i];

		// Model-based Marginal Variance
		double mb_varMarg  = 1 / mb_MV[i];
		double mb_statMarg = (betaMarg * betaMarg) / mb_varMarg;
		double mb_pvalMarg = (isnan(mb_statMarg) || mb_statMarg <= 0.0) ? NAN : boost::math::cdf(complement(chisq_dist_M, mb_statMarg));

		// Robust Marginal Variance
		double rb_varMarg  = 1 / rb_MV[i];
		double rb_statMarg = (betaMarg * betaMarg) / rb_varMarg;
		double rb_pvalMarg = (isnan(rb_statMarg) || rb_statMarg <= 0.0) ? NAN : boost::math::cdf(complement(chisq_dist_M, rb_statMarg));

		// Interaction effects
		subMatsubVecprod(mb_V, mb_U, betaInt, Sq1, Sq1, iss, is, 0);

		// Model-based Int P-value
		subMatrix(mb_V, VE, nExp, nExp, Sq1, nExp, iss + Sq1 + 1);
		matInv(VE, 1);
		matvecSprod(VE, betaInt, StempE, nExp, nExp, 1);
		double mb_statInt = 0.0;
		for (size_t j = 1; j < nExp1; j++) 
			mb_statInt += betaInt[j] * StempE[j-1];
		double mb_pvalInt = (isnan(mb_statInt) || mb_statInt <= 0.0) ? NAN : boost::math::cdf(complement(chisq_dist_Int, mb_statInt));

		// Model-based Joint P-value
		subMatrix(mb_V, Ai, nExp1, nExp1, nExp1, nExp1, iss);
		matInv(Ai, nExp1);
		matvecprod(Ai, betaInt, StempGE, nExp1, nExp1);
		double mb_statJoint = 0.0;
		for (size_t k = 0; k < nExp1; k++)
			mb_statJoint += betaInt[k] * StempGE[k];
		double mb_pvalJoint = (isnan(mb_statJoint) || mb_statJoint <= 0.0) ? NAN : boost::math::cdf(complement(chisq_dist_Joint, mb_statJoint));

		// Robust Int P-value
		subMatrix(rb_V, VE, nExp, nExp, Sq1, nExp, iss + Sq1 + 1);
		matInv(VE, 1);
		matvecSprod(VE, betaInt, StempE, nExp, nExp, 1);
		double rb_statInt = 0.0;
		for (size_t j = 1; j < nExp1; j++) 
			rb_statInt += betaInt[j] * StempE[j-1];
		double rb_pvalInt = (isnan(rb_statInt) || rb_statInt <= 0.0) ? NAN : boost::math::cdf(complement(chisq_dist_Int, rb_statInt));

		// Robust Joint P-value
		subMatrix(rb_V, Ai, nExp1, nExp1, nExp1, nExp1, iss);
		matInv(Ai, nExp1);
		matvecprod(Ai, betaInt, StempGE, nExp1, nExp1);
		double rb_statJoint = 0.0;
		for (size_t k = 0; k < nExp1; k++)
			rb_statJoint += betaInt[k] * StempGE[k];
		double rb_pvalJoint = (isnan(rb_statJoint) || rb_statJoint <= 0.0) ? NAN : boost::math::cdf(complement(chisq_dist_Joint, rb_statJoint));


		oss << snpid[i] << "\t" << chr[i] << "\t" << pos[i] << "\t" << nonEffectAllele[i] << "\t" << effectAllele[i] << "\t" << AF[i] << "\t" << nSamples[i] << "\t";
		oss << betaMarg << "\t" << mb_varMarg << "\t" << rb_varMarg << "\t";
		for (int j = 0; j < Sq1; j++)
		{
			oss << betaInt[j] << "\t";
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

		oss << mb_pvalMarg << "\t" << mb_pvalInt << "\t" << mb_pvalJoint << "\t";
		oss << rb_pvalMarg << "\t" << rb_pvalInt << "\t" << rb_pvalJoint << "\n";

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

	delete[] betaInt;
	delete[] StempE;
	delete[] StempGE;
	delete[] Ai;
	delete[] VE;
	delete[] mb_MU;
	delete[] rb_MU;
	delete[] mb_MV;
	delete[] rb_MV;
	delete[] mb_U;
	delete[] rb_U;
	delete[] mb_V;
	delete[] rb_V;
}
