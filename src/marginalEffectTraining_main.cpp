/*
Large-scale Genome Wide Association analysis (UKBiobank)
Copyright (C) 2018  Sheng Yang and Xiang Zhou

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <eigen3/Eigen/Dense>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <vector>
#include <cstring>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <map>
#include <set>
#include <algorithm>
#include <numeric>
#include <utility>
#include <sys/time.h>
#include <zlib.h>
#include "omp.h"
#include <math.h>
#include <random>

#include <boost/random.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/math/distributions/students_t.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/chi_squared.hpp>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include "Types.hpp"
#include "io.hpp"
#include "marginalTraining.hpp"

using namespace std;
using namespace Eigen;

int main(int argc, char *argv[]) {


	/*******************STEP1: Set parameters*******************/
	//Thresholds of MAF and INFO
	float bgenMinMAF = 1e-3;
	float bgenMinINFO = 0.8;
	float bgenMinHW = 1e-7;
	float bgenCallingRate = 0.99;

	//Number of threads 
	int threads = 50;

	//CV
	float prop = 0.8;
	const int seed = 20160522;

	cout << "STEP1 is ok!" << endl;

	/*******************STEP2: Load and process data*******************/
	//Read datasets
	char separate[] = ",";
	//Path of the R result: phenoStr, sqcNAStr and sqcStr
	string phenoStr = ;
	string sqcNAStr = ;
	string sqcStr = ;
	vector<vector <string> > pheno;
	vector<vector <string> > sqcNAIndex;
	vector<vector <string> > sqcData;

	readTable(phenoStr, 1, separate, pheno);
	readTable(sqcNAStr, 1, separate, sqcNAIndex);
	readTable(sqcStr, 1, separate, sqcData);
	cout << "All files are loaded!" << endl;

	//Process datasets
	vector <int> CVIndex;
	getCVIndex(prop, seed, sqcNAIndex, pheno, CVIndex);
		
	//Training data
	int NInitial = pheno.size();
	int sqcVar = sqcData[0].size() - 1;//delete the colname
	int NTraining = 0;
	vector<double> phenoTraining;
	vector<vector<double> > sqcTraining; 
	
	for (int i = 0; i < NInitial; i++) {

		if (CVIndex[i] == 1) {

			phenoTraining.push_back(atof(pheno[i][1].c_str()));
			vector <double> sqcTemp(sqcVar);
			for (int j = 0; j < sqcVar; j++) sqcTemp[j] = atof(sqcData[i][j + 1].c_str());
			sqcTraining.push_back(sqcTemp);
			NTraining++;
		}
	}
	
	VectorXd phenoTrainingVec = VectorXd::Map(&phenoTraining[0], NTraining);
	MatrixXd sqcTrainingMat(NTraining, sqcVar);
	for (int i = 0; i < NTraining; i++) sqcTrainingMat.row(i) = VectorXd::Map(&sqcTraining[i][0], sqcVar);

	cout << "STEP2 is ok!" << endl;

	/*******************STEP3: Get summary data and phypos data *******************/
	//Input files
	const int chrNum = 22;
	string chr[chrNum] = { "1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
						"11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22" };
	vector <string> bgen(chrNum);
	vector <string> output(chrNum);
	string startBgen = "/net/mulan/Biobank/rawdata/EGAD00010001225/001/ukb_imp_chr";
	string endBgen = "_v2.bgen";

	float threshold1 = 1e-3;
	float threshold2 = 1e-8;

	//Output files
	//Path of the results: outputSummary, outputInfoI and outputInfoII
	string outputSummary = "summaryData_cv_01.txt";
	string outputInfoI = "snpInfo(1e-3)_cv_01.txt";
	string outputInfoII = "snpInfo(1e-8)_cv_01.txt";
	ofstream foutSummary(outputSummary.c_str());
	ofstream foutInfoI(outputInfoI.c_str());
	ofstream foutInfoII(outputInfoII.c_str());

	for (int i = 0; i < chrNum; i++) {

		bgen[i] = startBgen + chr[i] + endBgen;

		cout << "Chr " << chr[i] << " is begining!" << endl;

		streamBgenTraining(chr[i], NTraining, bgen[i], phenoTrainingVec,
			sqcTrainingMat, CVIndex,
			bgenMinMAF, bgenMinINFO, bgenMinHW, bgenCallingRate,
			threads, threshold1, threshold2, foutSummary, foutInfoI, foutInfoII);

		cout << "Chr " << chr[i] << " is finished!" << endl;
	}

	foutSummary.close();
	foutInfoI.close();
	foutInfoII.close();

	cout << "STEP3 is ok!" << endl;
	return 0;
}
