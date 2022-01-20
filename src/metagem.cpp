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
	bool mb = cmd.mb;
	bool rb = cmd.rb;
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
		processFileHeader(Sq1, cmd.metaOpt, lc_intNames, cmd.fileNames[i], fip);
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
	std::vector<double> nSamples(nvars);
	std::vector<double> AF(nvars);
	

	// Get U vectors and V matrices from each file
	std::vector<double> mb_MU;
	std::vector<double> rb_MU;
	std::vector<double> mb_MV;
	std::vector<double> rb_MV;
	std::vector<double> mb_U;
	std::vector<double> rb_U;
	std::vector<double> mb_V;
	std::vector<double> rb_V;

	if (mb)
	{
		mb_MU.resize(uSize);
		mb_U.resize(uSize);
		mb_MV.resize(vSize);
		mb_V.resize(vSize);
	}

	if (rb)
	{
		rb_MU.resize(uSize);
		rb_U.resize(uSize);
		rb_MV.resize(vSize);
		rb_V.resize(vSize);
	}



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

		std::vector<int> betaIntColumn = fip->betaIntColumn[fileName];


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
		if (f == 0)
		{
			while(getline(file, line))
			{
				int ui = n * Sq1;
				int vi = n * Sq1 * Sq1;

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

				std::transform(a1.begin(), a1.end(), a1.begin(), [](char c){ return std::tolower(c); });
				std::transform(a2.begin(), a2.end(), a2.begin(), [](char c){ return std::tolower(c); });
				nonEffectAllele[n] = a1;
				effectAllele[n]    = a2;



				// Beta Int
				std::vector<double> fileBetaInt(Sq1);
				for (int i = 0; i < Sq1; i++) {
					fileBetaInt[i] = std::stod(values[betaIntColumn[i]]);
				}

				// Model-based
				if (mb)
				{		
					int mb_seMargColumn = fip->mb_seMargColumn[fileName];
					std::vector<int> mb_covIntColumn = fip->mb_covIntColumn[fileName];
					
					double mb_fileBM = std::stod(values[betaMargColumn]);
					double mb_fileMV = std::stod(values[mb_seMargColumn]);

					for (int i = 0; i < Sq1; i++) {
						int isq = (i * Sq1);
						int ii  = vi + (i * Sq1);
						for (int j = 0; j < Sq1; j++) {
							mb_V[ii + j] = std::stod(values[mb_covIntColumn[isq + j]]);
						}
						mb_V[ii + i] *= mb_V[ii + i];
					}

					mb_MV[n] += (1 / (mb_fileMV * mb_fileMV));
					mb_MU[n] += (mb_MV[n] * mb_fileBM);

					subMatInv(&mb_V[0], Sq1, vi);
					subMatsubVecprod(&mb_V[0], &fileBetaInt[0], &mb_U[0], Sq1, Sq1, vi, 0, ui);
				}

				// Robust
				if (rb)
				{
					int rb_seMargColumn = fip->rb_seMargColumn[fileName];
					std::vector<int> rb_covIntColumn = fip->rb_covIntColumn[fileName];

					double rb_fileBM = std::stod(values[betaMargColumn]);
					double rb_fileMV = std::stod(values[rb_seMargColumn]);

					for (int i = 0; i < Sq1; i++) {
						int isq = (i * Sq1);
						int ii  = vi + (i * Sq1);
						for (int j = 0; j < Sq1; j++) {
							rb_V[ii + j] = std::stod(values[rb_covIntColumn[isq + j]]);
						}
						rb_V[ii + i] *= rb_V[ii + i];
					}
	
					rb_MV[n] += (1 / (rb_fileMV * rb_fileMV));
					rb_MU[n] += (rb_MV[n] * rb_fileBM);
					subMatInv(&rb_V[0], Sq1, vi);
					subMatsubVecprod(&rb_V[0], &fileBetaInt[0], &rb_U[0], Sq1, Sq1, vi, 0, ui);
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

				std::string snpName = values[snpColumn];

				std::vector<double> mb_fileU(Sq1);
				std::vector<double> mb_fileV(Sq1 * Sq1);
				std::vector<double> rb_fileU(Sq1);
				std::vector<double> rb_fileV(Sq1 * Sq1);
				if (mb)
				{
					std::vector<int> mb_covIntColumn = fip->mb_covIntColumn[fileName];

					for (int i = 0; i < Sq1; i++) {
						int ii = i * Sq1;
						for (int j = 0; j < Sq1; j++) {
							mb_fileV[ii + j] = std::stod(values[mb_covIntColumn[ii + j]]);
						}
						mb_fileV[ii + i] *= mb_fileV[ii + i];
					}
				}

				if (rb)
				{
					std::vector<int> rb_covIntColumn = fip->rb_covIntColumn[fileName];

					for (int i = 0; i < Sq1; i++) {
						int ii = i * Sq1;
						for (int j = 0; j < Sq1; j++) {
							rb_fileV[ii + j] = std::stod(values[rb_covIntColumn[ii + j]]);
						}
						rb_fileV[ii + i] *= rb_fileV[ii + i];
					}
				}


				std::unordered_map<std::string, int>::iterator it = snpid_idx.find(snpName);
				if (it != snpid_idx.end())
				{
					int val = it->second;

					double dir = 1.0;
					int flip   = flipAllele(nonEffectAllele[val], effectAllele[val], values[nonEffectColumn], values[effectColumn]);

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
						n++;
						continue;
					}

					// Beta Int
					std::vector<double> fileBetaInt(Sq1);
					for (int i = 0; i < Sq1; i++) {
						fileBetaInt[i] = std::stod(values[betaIntColumn[i]]) * dir;
					}

					int ui = val * Sq1;
					int vi = val * (Sq1 * Sq1);

					if (mb)
					{
						int mb_seMargColumn = fip->mb_seMargColumn[fileName];
						double mb_fileBM = std::stod(values[betaMargColumn]);
						double mb_fileMV = std::stod(values[mb_seMargColumn]);
						
						mb_fileMV = 1 / (mb_fileMV * mb_fileMV);
						mb_MV[val] += mb_fileMV;
						mb_MU[val] += (mb_fileMV * (mb_fileBM * dir));

						matInv(&mb_fileV[0], Sq1);
						matvecprod(&mb_fileV[0], &fileBetaInt[0], &mb_fileU[0], Sq1, Sq1);
						for (int i = 0; i < Sq1; i++)
						{
							int ii = i * Sq1;
							mb_U[ui+ i] += mb_fileU[i];
							for (int j = 0; j < Sq1; j++)
							{
								mb_V[vi + ii + j] += mb_fileV[ii + j];
							}
						}
					}

					if (rb)
					{
						int rb_seMargColumn = fip->rb_seMargColumn[fileName];
						double rb_fileBM = std::stod(values[betaMargColumn]);
						double rb_fileMV = std::stod(values[rb_seMargColumn]);

						rb_fileMV = 1 / (rb_fileMV * rb_fileMV);
						rb_MV[val] += rb_fileMV;
						rb_MU[val] += (rb_fileMV * (rb_fileBM * dir));

						matInv(&rb_fileV[0], Sq1);
						matvecprod(&rb_fileV[0], &fileBetaInt[0], &rb_fileU[0], Sq1, Sq1);
						for (int i = 0; i < Sq1; i++)
						{
							int ii = i * Sq1;
							rb_U[ui+ i] += rb_fileU[i];
							for (int j = 0; j < Sq1; j++)
							{
								rb_V[vi + ii + j] += rb_fileV[ii + j];
							}
						}
					}

					n++;
				}
				else {
					snpid.push_back(snpName);
					chr.push_back(values[chrColumn]);
					pos.push_back(values[posColumn]);

					std::string a1 = values[nonEffectColumn];
					std::string a2 = values[effectColumn];

					std::transform(a1.begin(), a1.end(), a1.begin(), [](char c){ return std::tolower(c); });
					std::transform(a2.begin(), a2.end(), a2.begin(), [](char c){ return std::tolower(c); });
					nonEffectAllele.push_back(a1);
					effectAllele.push_back(a2);

					double ns = std::stod(values[nSampleColumn]);
					nSamples.push_back(ns);
 					AF.push_back(std::stod(values[afColumn]) * 2.0 * ns);

					// Beta Int
					std::vector<double> fileBetaInt(Sq1);
					for (int i = 0; i < Sq1; i++) {
						fileBetaInt[i] = std::stod(values[betaIntColumn[i]]);
					}

					if (mb)
					{
						int mb_seMargColumn = fip->mb_seMargColumn[fileName];
						double mb_fileBM = std::stod(values[betaMargColumn]);
						double mb_fileMV = std::stod(values[mb_seMargColumn]);

						mb_fileMV = 1 / (mb_fileMV * mb_fileMV);
						mb_MV.push_back(mb_fileMV);
						mb_MU.push_back(mb_fileMV * mb_fileBM);

						matInv(&mb_fileV[0], Sq1);
						matvecprod(&mb_fileV[0], &fileBetaInt[0], &mb_fileU[0], Sq1, Sq1);
						mb_V.insert(mb_V.end(), mb_fileV.begin(), mb_fileV.end());
						mb_U.insert(mb_U.end(), mb_fileU.begin(), mb_fileU.end());
					}

					if (rb)
					{
						int rb_seMargColumn = fip->rb_seMargColumn[fileName];
						double rb_fileBM = std::stod(values[betaMargColumn]);
						double rb_fileMV = std::stod(values[rb_seMargColumn]);

						rb_fileMV = 1 / (rb_fileMV * rb_fileMV);
						rb_MV.push_back(rb_fileMV);
						rb_MU.push_back(rb_fileMV * rb_fileBM);

						matInv(&rb_fileV[0], Sq1);
						matvecprod(&rb_fileV[0], &fileBetaInt[0], &rb_fileU[0], Sq1, Sq1);
						rb_V.insert(rb_V.end(), rb_fileV.begin(), rb_fileV.end());
						rb_U.insert(rb_U.end(), rb_fileU.begin(), rb_fileU.end());
					}

					snpid_idx[snpName] = nvars;
					nvars++;
					n++;
				}
			}

		}
		file.close();

		if (n != nrow) {
			cerr << "\nERROR: The number of variants read in [" << fileName << "] does not match original count.\n\n";
			exit(1);
		}
	}

	// Create the output file and write column header names
	printOutputHeader(cmd.outFile, cmd.metaOpt, Sq1, cmd.intNames);
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

	double* mb_betaInt = new double[Sq1];
	double* rb_betaInt = new double[Sq1];
	double* StempE  = new double[nExp];
	double* StempGE = new double[nExp1];
	double* Ai      = new double[nExp1 * nExp1];
	double* VE      = new double[nExp * nExp];


	for (int i = 0; i < nvars; i++)
	{

		int is  = (i * Sq1);
		int iss = (i * Sq1 * Sq1);

		AF[i] = AF[i] / 2.0 / nSamples[i];

		if (mb)
		{
			subMatInv(&mb_V[0], Sq1, iss);

			// Model-based Marginal Variance
			mb_betaMarg = mb_MU[i] / mb_MV[i];
			mb_varMarg  = 1 / mb_MV[i];
			double mb_statMarg = (mb_betaMarg * mb_betaMarg) / mb_varMarg;
			mb_pvalMarg = (std::isnan(mb_statMarg) || mb_statMarg <= 0.0) ? NAN : boost::math::cdf(complement(chisq_dist_M, mb_statMarg));

			// Interaction effects
			subMatsubVecprod(&mb_V[0], &mb_U[0], mb_betaInt, Sq1, Sq1, iss, is, 0);

			// Model-based Int P-value
			subMatrix(&mb_V[0], VE, nExp, nExp, Sq1, nExp, iss + Sq1 + 1);
			matInv(VE, 1);
			matvecSprod(VE, mb_betaInt, StempE, nExp, nExp, 1);
			double mb_statInt = 0.0;
			for (size_t j = 1; j < nExp1; j++) 
				mb_statInt += mb_betaInt[j] * StempE[j-1];
			mb_pvalInt = (std::isnan(mb_statInt) || mb_statInt <= 0.0) ? NAN : boost::math::cdf(complement(chisq_dist_Int, mb_statInt));
	
			// Model-based Joint P-value
			subMatrix(&mb_V[0], Ai, nExp1, nExp1, nExp1, nExp1, iss);
			matInv(Ai, nExp1);
			matvecprod(Ai, mb_betaInt, StempGE, nExp1, nExp1);
			double mb_statJoint = 0.0;
			for (size_t k = 0; k < nExp1; k++)
				mb_statJoint += mb_betaInt[k] * StempGE[k];
			mb_pvalJoint = (std::isnan(mb_statJoint) || mb_statJoint <= 0.0) ? NAN : boost::math::cdf(complement(chisq_dist_Joint, mb_statJoint));
	
		}

		if (rb)
		{
			subMatInv(&rb_V[0], Sq1, iss);

			// Robust Marginal Variance
			rb_betaMarg = rb_MU[i] / rb_MV[i];
			rb_varMarg  = 1 / rb_MV[i];
			double rb_statMarg = (rb_betaMarg * rb_betaMarg) / rb_varMarg;
			rb_pvalMarg = (std::isnan(rb_statMarg) || rb_statMarg <= 0.0) ? NAN : boost::math::cdf(complement(chisq_dist_M, rb_statMarg));

			// Interaction effects
			subMatsubVecprod(&rb_V[0], &rb_U[0], rb_betaInt, Sq1, Sq1, iss, is, 0);

			// Robust Int P-value
			subMatrix(&rb_V[0], VE, nExp, nExp, Sq1, nExp, iss + Sq1 + 1);
			matInv(VE, 1);
			matvecSprod(VE, rb_betaInt, StempE, nExp, nExp, 1);
			double rb_statInt = 0.0;
			for (size_t j = 1; j < nExp1; j++) 
				rb_statInt += rb_betaInt[j] * StempE[j-1];
			rb_pvalInt = (std::isnan(rb_statInt) || rb_statInt <= 0.0) ? NAN : boost::math::cdf(complement(chisq_dist_Int, rb_statInt));

			// Robust Joint P-value
			subMatrix(&rb_V[0], Ai, nExp1, nExp1, nExp1, nExp1, iss);
			matInv(Ai, nExp1);
			matvecprod(Ai, rb_betaInt, StempGE, nExp1, nExp1);
			double rb_statJoint = 0.0;
			for (size_t k = 0; k < nExp1; k++)
				rb_statJoint += rb_betaInt[k] * StempGE[k];
			rb_pvalJoint = (std::isnan(rb_statJoint) || rb_statJoint <= 0.0) ? NAN : boost::math::cdf(complement(chisq_dist_Joint, rb_statJoint));
		}


		oss << snpid[i] << "\t" << chr[i] << "\t" << pos[i] << "\t" << nonEffectAllele[i] << "\t" << effectAllele[i] << "\t" << nSamples[i] << "\t" << AF[i] << "\t";

		// Marginal summary stats
		if (mb)
		{
			oss << mb_betaMarg << "\t" << mb_varMarg << "\t";
		}
		if (rb)
		{
			oss << rb_betaMarg << "\t" << rb_varMarg << "\t";
		}

		// Model based summary stats
		if (mb) {
			for (int j = 0; j < Sq1; j++)
			{
				oss << mb_betaInt[j] << "\t";
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
		}

		// Robust summary stats
		if (rb) 
		{
			for (int j = 0; j < Sq1; j++)
			{
				oss << rb_betaInt[j] << "\t";
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
		}

		if (mb)
		{
			oss << mb_pvalMarg << "\t" << mb_pvalInt << "\t" << mb_pvalJoint << ((rb) ? "\t" : "\n");
		}

		if (rb)
		{
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

	cout << "Done." << endl;
}
