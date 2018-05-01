#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <boost/random.hpp>

#include "io.hpp"
#include "Types.hpp"

using namespace std;

/*************** Load files (csv, txt) ***************/
//csv: separator = ","
//txt: separator = "\t"
int readTable(string infile, int skip, char *separator,
			vector< vector<string> > &data)
{

	string row, element;
	ifstream fileStream(infile.c_str());

	if (!fileStream) {

		cerr << "ERROR: " << infile << " dose not exist!" << endl;
	}

	int i = 0;
	while (getline(fileStream, row)) {

		vector<string> rowVec;
		stringstream rowStream(row);
		if (i > skip - 1) {

			while (getline(rowStream, element, *separator)) rowVec.push_back(element);
			data.push_back(rowVec);
		}

		i++;
	}

	fileStream.close();
	return 0;
}

/*************** Transform base to code ***************/
int baseToCode(string base)
{

	int code;
	if (base == "A") {
		code = 0;
	}
	else if (base == "T") {
		code = 1;
	}
	else if (base == "C") {
		code = 2;
	}
	else if (base == "G") {
		code = 3;
	}
	else {
		code = 4;
	}
	return code;
}

/*************** Transform code to base ***************/
//output the ref allele and alt allele
string codeToBase(int code)
{

	string base;
	if (code == 0) {
		base = 'A';
	}
	else if (code == 1) {
		base = 'T';
	}
	else if (code == 2) {
		base = 'C';
	}
	else if (code == 3) {
		base = 'G';
	}
	else if (code == 4) {
		base = 'X';
	}
	return base;
}

/*************** output the summary data as the format of ldpred ***************/
string standardOutputSummaryData(vector <double> snpInfo, string snpName)
{

	ostringstream fout;
	int chrNum = (int)snpInfo[0];
	string chr = "chr" + to_string(chrNum);
	int32 pos = (int32)snpInfo[1];
	double reffreq = 1 - snpInfo[4];
	string ref, alt;
	ref = codeToBase(snpInfo[2]);
	alt = codeToBase(snpInfo[3]);

	fout << chr << "\t" << to_string(pos) << "\t" << ref << "\t" << alt << "\t" << reffreq << "\t" << snpInfo[5] << "\t" << snpName << "\t" << snpInfo[6] << "\t" << snpInfo[7];
	fout << endl;

	return fout.str();
}

/*************** output the snp information ***************/
//chr, phypos, pval, filepos
string standardOutputInfo(vector <double> snpInfo)
{

	ostringstream fout;
	int chrNum = (int)snpInfo[0];
	string chr = to_string(chrNum);
	int32 phypos = (int32)snpInfo[1];
	int64 filepos = (int64)snpInfo[3];

	fout << chr << "\t" << to_string(phypos) << "\t" << snpInfo[2] << "\t" << filepos;
	fout << endl;

	return fout.str();
}

/*************** get the index of training samples ***************/
int getCVIndex(float prop, int seed, vector<vector <string> > sqcNAIndex, vector<vector <string> > pheno,
						vector <int> &CVIndex)
{

	int initialN = pheno.size();
	cout << "Sample size (genotype success): " << initialN << endl;

	//sqc NA
	int sqcN = 0;
	for (int i = 0; i < sqcNAIndex.size(); i++) {

		if (sqcNAIndex[i][0] == "1") sqcN++;
	}
	cout << "Sample size (sqc success): " << sqcN << endl;

	//Pheontype NA
	vector <int> phenoNAIndex(initialN);
	fill(phenoNAIndex.begin(), phenoNAIndex.end(), 1);

	for (int i = 0; i < initialN; i++) {

		if (pheno[i][1] == "NA") {

			phenoNAIndex[i] = 0;
		}
	}
	int NA_N = accumulate(phenoNAIndex.begin(), phenoNAIndex.end(), 0);
	cout << "Sample size (not NA phenotype): " << NA_N << endl;

	//Summarize sqc NA and phenotype NA
	vector <int> NAIndex(initialN);
	fill(NAIndex.begin(), NAIndex.end(), 0);
	for (int i = 0; i < initialN; i++) {

		if (phenoNAIndex[i] == 1 && sqcNAIndex[i][0] == "1")

			NAIndex[i] = 1;
	}

	const int nlyN = accumulate(NAIndex.begin(), NAIndex.end(), 0);

	cout << "Sample size (analysis): " << nlyN << endl;

	//Cross validation
	default_random_engine generator(seed);
	binomial_distribution<int> binDistribution(1, prop);
	for (int i = 0; i < initialN; i++) {

		if (NAIndex[i] == 1) {

			CVIndex.push_back(binDistribution(generator));
		}
		else {

			CVIndex.push_back(-9);
		}
	}

	int trainingN = 0;
	for (int i = 0; i < initialN; i++) {

		if (CVIndex[i] == 1) {
			trainingN++;
		}
	}
	cout << "Sample size (training): " << trainingN << endl;

	return 0;
}