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
#include <vector>
#include <string>
#include <iostream>

#include "lmFunc.hpp"
#include "snpEff.hpp"
#include "mathFunc.hpp"
#include "steprd.hpp"

using namespace std;
using namespace Eigen;

STEPRD::STEPRD(void) :
	version("1.0"), date("05/17/2018"), year("2018")
{}

void STEPRD::printHeader(void)
{
	cout << endl;
	cout << "*************************************************************" << endl;
	cout << "  Large-scale Genome Wide Association analysis (UKBiobank)   " << endl;
	cout << "  Version " << version << ", " << date << "                  " << endl;
	cout << "  Visit http://www.xzlab.org/software.html For Updates       " << endl;
	cout << "  (C) " << year << " Sheng Yang, Xiang Zhou                  " << endl;
	cout << "  GNU General Public License                                 " << endl;
	cout << "  For Help, Type ./steprd -h                                 " << endl;
	cout << "*************************************************************" << endl;
	cout << endl;

	return;
}

void STEPRD::printHelp(void) {

	cout << " FILE I/O RELATED OPTIONS" << endl;
	cout << " -prop      [num]       " << " specify input the proportion of training data." << endl;
	cout << " -seed      [num]       " << " specify input the seed of cross validation." << endl;
	cout << " -t1        [num]       " << " specify input the threashold 1 (1e-8)." << endl;
	cout << " -t2        [num]       " << " specify input the threashold 2 (1e-3)." << endl;
	cout << " -ts        [num]       " << " specify input the threashold of stepwise selection." << endl;
	cout << " -pheno     [filename]  " << " specify input the phenotype data." << endl;
	cout << " -sqc       [filename]  " << " specify input the sqc NA label." << endl;
	cout << " -summ      [filename]  " << " specify input the summary data." << endl;
	cout << "             includes: chr, position, P value and position of file" << endl;
	cout << " -len1      [num]       " << " specify input the distance between snps less than t1." << endl;
	cout << " -len2      [num]       " << " specify input the distance between snps less than t2." << endl;
	cout << " -thread    [num]       " << " specify input the applied threads." << endl;
	cout << " -bgen      [filename]  " << " specify input the bgen data." << endl;
	cout << "             requires: *.bgen files" << endl;
	cout << " -outpath   [path]      " << " specify input the output path." << endl;
	cout << " -snpfile   [filename]  " << " specify input the filename of summary data." << endl;
	cout << " -phenofile [filename]  " << " specify input the filename of testing and prediction phenotype data." << endl;

	return;
}

void STEPRD::Assign(int argc, char ** argv, PARAM &cPar) {

	string str;
	for (int i = 0; i < argc; i++) {

		if (strcmp(argv[i], "--PROP") == 0 || strcmp(argv[i], "-prop") == 0) {

			if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.prop = atof(str.c_str());
		}
		else if (strcmp(argv[i], "--SEED") == 0 || strcmp(argv[i], "-seed") == 0) {

			if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.seed = atoi(str.c_str());
		}
		else if (strcmp(argv[i], "--PHENO") == 0 || strcmp(argv[i], "-pheno") == 0) {

			if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.phenoStr = str;
		}
		else if (strcmp(argv[i], "--SQC") == 0 || strcmp(argv[i], "-sqc") == 0) {

			if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.sqcNAStr = str;
		}
		else if (strcmp(argv[i], "--SUMM") == 0 || strcmp(argv[i], "-summ") == 0) {

			if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.summDataStr = str;
		}
		else if (strcmp(argv[i], "--T1") == 0 || strcmp(argv[i], "-t1") == 0) {

			if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.threshold1 = atof(str.c_str());
		}
		else if (strcmp(argv[i], "--T2") == 0 || strcmp(argv[i], "-t2") == 0) {

			if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.threshold2 = atof(str.c_str());
		}
		else if (strcmp(argv[i], "--TS") == 0 || strcmp(argv[i], "-ts") == 0) {

			if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.thresholdStep = atof(str.c_str());
		}
		else if (strcmp(argv[i], "--LEN1") == 0 || strcmp(argv[i], "-len1") == 0) {

			if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.Len1 = atof(str.c_str());
		}
		else if (strcmp(argv[i], "--LEN2") == 0 || strcmp(argv[i], "-len2") == 0) {

			if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.Len2 = atof(str.c_str());
		}
		else if (strcmp(argv[i], "--THREAD") == 0 || strcmp(argv[i], "-thread") == 0) {

			if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.threads = atoi(str.c_str());
		}
		else if (strcmp(argv[i], "--BGENFILE") == 0 || strcmp(argv[i], "-bgenfile") == 0) {

			if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.bgenStart = str;
		}
		else if (strcmp(argv[i], "--OUTPATH") == 0 || strcmp(argv[i], "-outpath") == 0) {

			if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.outpath = str;
		}
		else if (strcmp(argv[i], "--SNPFILE") == 0 || strcmp(argv[i], "-snpfile") == 0) {

			if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.snpfile = str;
		}
		else if (strcmp(argv[i], "--PHENOFILE") == 0 || strcmp(argv[i], "-phenofile") == 0) {

			if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.phenofile = str;
		}
	}
	return;
}

void STEPRD::BatchRun(PARAM &cPar) {

	  ////////////////////////////////////////////
	 ///         Buile stepwise model         ///
	////////////////////////////////////////////

	STEPSEL SS;
	//Load data
	cout << "Summary data input ... " << endl;
	char separate1[] = ",";
	char separate2[] = "\t";
	vector<vector <string> > pheno;
	vector<vector <string> > sqcNAIndex;
	vector<vector <string> > summData;
	readTable(cPar.phenoStr, 1, separate1, pheno);
	readTable(cPar.sqcNAStr, 1, separate1, sqcNAIndex);
	readTable(cPar.summDataStr, 1, separate2, summData);

	//Process phenotype datasets: get the training and test phenotype
	vector <int> CVIndex;
	getCVIndex(cPar.prop, cPar.seed, sqcNAIndex, pheno, CVIndex);

	//Process phenotype data
	cout << "Process phenotype data ..." << endl;
	int NInitial = pheno.size();
	vector<double> phenoTraining;
	vector<double> phenoTesting;
	for (int i = 0; i < NInitial; i++) {

		if (CVIndex[i] == 1) {

			phenoTraining.push_back(atof(pheno[i][0].c_str()));
		}
		if (CVIndex[i] == 0) {

			phenoTesting.push_back(atof(pheno[i][0].c_str()));
		}
	}
	int NTraining = phenoTraining.size();
	int NTesting = phenoTesting.size();
	VectorXd phenoTrainingVec = VectorXd::Map(&phenoTraining[0], NTraining);
	double phenoTrainingMean = phenoTrainingVec.mean();
	double phenoTrainingSd = calcSd(phenoTrainingVec);
	VectorXd phenoTestingStd = standardizeVec(VectorXd::Map(&phenoTesting[0], NTesting), phenoTrainingMean, phenoTrainingSd);

	//Get the training genotype of significant SNPs (t1) in 2Mb
	cout << "Get independent SNPs (<T1) of training data ..." << endl;
	string bgenEnd = "_v2.bgen";
	MatrixXd genoTrainingStdTh1;
	vector <vector <double> > SummTh1Sig;
	SS.getGenotype1(summData, cPar.Len1, 0, cPar.threshold1, NTraining, NTesting, CVIndex, cPar.bgenStart, bgenEnd,
		cPar.threads, "Training", genoTrainingStdTh1, SummTh1Sig);
	SummTh1Sig.resize(0);
	// cout << SummTh1Sig.size() << endl;
	VectorXd resid = mvlmResid(genoTrainingStdTh1, phenoTrainingStd);
	//for (int i = 0; i < SummTh1Sig.size(); i++) {
	//	cout << SummTh1Sig[i][0] << "\t" << SummTh1Sig[i][1] << "\t" << SummTh1Sig[i][2] << "\t" << SummTh1Sig[i][3] << "\t" << endl;
	//}

	// Calculate P value of genoTrainingTh2
	cout << "Calculate and sort P value of independent SNPs (<T2) of training data ..." << endl;
	vector <double> pTot;
	vector <vector <double> > snpIdp;
	SS.getGenotype2(summData, cPar.Len2, cPar.threshold1, cPar.threshold2, NTraining, NTesting, CVIndex, cPar.bgenStart, bgenEnd, 
					cPar.threads, "Training", resid, pTot, snpIdp);

	vector <vector <double> > pTotIdx(pTot.size(), vector<double>(2));
	for (int i = 0; i < pTot.size(); i++) {

		pTotIdx[i][0] = i;
		pTotIdx[i][1] = pTot[i];
		//cout << pTotIdx[i][0] << "\t" << pTotIdx[i][1] << endl;
	}
	sort(pTotIdx.begin(), pTotIdx.end(), SortCol1);
	// for (int i = 0; i < 200; i++) cout << pTotIdx[i][0] << "\t" << pTotIdx[i][1] << endl;
	
	//Build stepwise model
	cout << "Build stepwise model ... " << endl;
	vector <vector <double> > SummTh2Sig;
	MatrixXd genoTrainingStd;
	lmStepwise(genoTrainingStdTh1, phenoTrainingStd, snpIdp, CVIndex, cPar.thresholdStep, cPar.bgenStart, bgenEnd, pTotIdx, SummTh2Sig, genoTrainingStd);
	VectorXd betaAll = genoTrainingStd.transpose() * phenoTrainingStd / NTraining;
	VectorXd phenoTrainingHat = genoTrainingStd * betaAll;
	double trainingMSE = calcMSE(phenoTrainingStd, phenoTrainingHat);
	cout << "MSE for training data is: " << trainingMSE << endl;
	cout << "Stepwise model is ok!" << endl;
	genoTrainingStdTh1.resize(0, 0);
	genoTrainingStd.resize(0, 0);

	//Get independent SNPs (<t1) of testing data
	cout << "Get independent SNPs (<T1) of testing data ..." << endl;
	MatrixXd genoTestingStdTh1;
	// vector <vector <double> > SummTh1Sig;
	SS.getGenotype1(summData, cPar.Len1, 0, cPar.threshold1, NTraining, NTesting, CVIndex, cPar.bgenStart, bgenEnd, cPar.threads, "Testing", 
		genoTestingStdTh1, SummTh1Sig);

	//Calculate P value of genoTrainingTh2
	cout << "Signigifcant genoTestingTh2 ..." << endl;
	vector <vector <double> > genoTestingTh2;
	for (int i = 0; i < SummTh2Sig.size(); i++) {

		vector <double> genoTestingVec;
		SS.streamBgenSingle(SummTh2Sig[i], CVIndex, cPar.bgenStart, bgenEnd, "Testing", genoTestingVec);
		genoTestingTh2.push_back(genoTestingVec);
	}
	//transformate vector to MatrixXd
	int varNumTh2 = genoTestingTh2.size();
	MatrixXd genoTestingStdTh2(varNumTh2, NTesting);
	for (int i = 0; i < varNumTh2; i++) {

		genoTestingStdTh2.row(i) = standardizeVec(VectorXd::Map(&genoTestingTh2[i][0], NTesting), SummTh2Sig[i][4], SummTh2Sig[i][5]);
	}
	genoTestingStdTh2.transposeInPlace();
	// cout << genoTestingStdTh2.cols() << "\t" << genoTestingStdTh2.rows() << endl;

	//Fit model of testing data
	cout << "Fit model of testing data ... " << endl;
	MatrixXd genoTestingMatComb(NTesting, genoTestingStdTh1.cols() + genoTestingStdTh2.cols());
	genoTestingMatComb << genoTestingStdTh1, genoTestingStdTh2;
	VectorXd phenoHat = genoTestingMatComb * betaAll;
	double MSE = calcMSE(phenoTestingStd, phenoHat);
	double COR = calcCOR(phenoTestingStd, phenoHat);
	cout << "STEPRD MSE for the whole genome is: " << MSE << endl;
	cout << "STEPRD COR for the whole genome is: " << COR << endl;

	//Integrate SNP information
	vector <vector <double> > SummSig;
	SummSig = SummTh1Sig;
	for (int i = 0; i < SummTh2Sig.size(); i++) {

		SummSig.push_back(SummTh2Sig[i]);
	}
	//cout << SummTh2Sig.size() << endl;
	//cout << SummSig.size() << endl;
	//cout << betaAll.size() << endl;
	
	//Output SNP information and betaHat
	string snpOutput = cPar.outpath + "/" + cPar.snpfile;
	ofstream snpFout(snpOutput.c_str());
	snpFout << "CHR" << "\t" << "POS" << "\t" << "P" << "\t" << "FPOS" << "\t" << "beta" << endl;
	for (int i = 0; i < SummSig.size(); i++) {

		snpFout << (int)SummSig[i][0] << "\t" << (uint)SummSig[i][1] << "\t" << SummSig[i][2] << "\t" << (uint)SummSig[i][3] << "\t" << betaAll(i);
		snpFout << endl;
	}
	
	//Output pehnotype 
	string phenoOutput = cPar.outpath + "/" + cPar.phenofile;
	ofstream phenoFout(phenoOutput.c_str());
	phenoFout << "y" << "\t" << "yHat" << endl;
	for (int i = 0; i < SummSig.size(); i++) {
	
		phenoFout << phenoTestingStd(i) << "\t" << phenoHat(i) << endl;
	}
	return;
}