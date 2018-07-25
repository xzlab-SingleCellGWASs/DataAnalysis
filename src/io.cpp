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

		if (pheno[i][0] == "NA") {

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
	//nlyN = accumulate(NAIndex.begin(), NAIndex.end(), 0);
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

	//CVIndex: -9, 0, 1
	int trainingN = 0;
	for (int i = 0; i < initialN; i++) {

		if (CVIndex[i] == 1) {
			trainingN++;
		}
	}
	cout << "Sample size (training): " << trainingN << endl;

	return 0;
}