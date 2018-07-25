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

#include <vector>
#include <string>
#include <iostream>

#include "Types.hpp"
#include "io.hpp"
#include "snpEff.hpp"
#include "mathFunc.hpp"
#include "summ.hpp"

using namespace std;

SUMM::SUMM(void) :
version("1.0"), date("05/17/2018"), year("2018")
{}

void SUMM::printHeader(void)
{
	cout << endl;
	cout << "************************************************************" << endl;
	cout << "  Large-scale Genome Wide Association analysis (UKBiobank)  " << endl;
	cout << "  Version " << version << ", " << date << "                 " << endl;
	cout << "  Visit http://www.xzlab.org/software.html For Updates      " << endl;
	cout << "  (C) " << year << " Sheng Yang, Xiang Zhou                 " << endl;
	cout << "  GNU General Public License                                " << endl;
	cout << "  For Help, Type ./summ -h                                  " << endl;
	cout << "************************************************************" << endl;
	cout << endl;

	return;
}

void SUMM::printHelp(void) {

	cout << " FILE I/O RELATED OPTIONS" << endl;
	cout << " -maf       [num]       " << " specify input the minimum of maf." << endl;
	cout << " -info      [num]       " << " specify input the minimum of snp infomation." << endl;
	cout << " -hwe       [num]       " << " specify input the minimum of p value of Hardy-Weinberg Equilibrium." << endl;
	cout << " -call      [num]       " << " specify input the minimum of snp calling rate." << endl;
	cout << " -thread    [num]       " << " specify input the applied threads." << endl;
	cout << " -prop      [num]       " << " specify input the proportion of training data." << endl;
	cout << " -seed      [num]       " << " specify input the seed of cross validation." << endl;
	cout << " -pheno     [filename]  " << " specify input the phenotype data." << endl;
	cout << " -sqc       [filename]  " << " specify input the sqc NA label." << endl;
	cout << " -bgen      [filename]  " << " specify input the bgen data." << endl;
	cout << "             requires: *.bgen files" << endl;
	cout << " -outpath   [path]      " << " specify input the output path." << endl;
	cout << " -outfile   [filename]  " << " specify input the filename of summary data." << endl;
	cout << " -chr       [num]       " << " specify input the chromosome." << endl;
	cout << " -cv        [num]       " << " specify input the time of cross validation." << endl;
	return;
}

void SUMM::Assign(int argc, char ** argv, PARAM &cPar) {

	string str;
	for (int i = 0; i < argc; i++) {

		if (strcmp(argv[i], "--MAF") == 0 || strcmp(argv[i], "-maf") == 0) {

			if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.bgenMinMAF = atof(str.c_str());
		}
		else if (strcmp(argv[i], "--INFO") == 0 || strcmp(argv[i], "-info") == 0) {
			
			if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.bgenMinINFO = atof(str.c_str());
		}
		else if (strcmp(argv[i], "--HWE") == 0 || strcmp(argv[i], "-hwe") == 0) {
			
			if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.bgenMinHWE = atof(str.c_str());
		}
		else if (strcmp(argv[i], "--CALL") == 0 || strcmp(argv[i], "-call") == 0) {
			
			if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.bgenMinCallingRate = atof(str.c_str());
		}
		else if (strcmp(argv[i], "--THREAD") == 0 || strcmp(argv[i], "-thread") == 0) {

			if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.threads = atof(str.c_str());
		}
		else if (strcmp(argv[i], "--PROP") == 0 || strcmp(argv[i], "-prop") == 0) {

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
		else if (strcmp(argv[i], "--BGEN") == 0 || strcmp(argv[i], "-bgen") == 0) {

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
		else if (strcmp(argv[i], "--OUTFILE") == 0 || strcmp(argv[i], "-outfile") == 0) {

			if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.outfile = str;
		}
		else if (strcmp(argv[i], "--CHR") == 0 || strcmp(argv[i], "-chr") == 0) {

			if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.chr = str;
		}
		else if (strcmp(argv[i], "--CV") == 0 || strcmp(argv[i], "-cv") == 0) {

			if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.cv = str;
		}
	}
	return;
}

void SUMM::BatchRun(PARAM &cPar) {

	MKSUMM MS;

	//Read Files
	cout << "sqc Data ... " << endl;
	char separate[] = ",";
	vector<vector <string> > sqcNAIndex;
	readTable(cPar.sqcNAStr, 1, separate, sqcNAIndex);

	//Phenotype Data
	cout << "Phenotype Data ... " << endl;
	vector<vector <string> > pheno;
	readTable(cPar.phenoStr, 1, separate, pheno);

	//Process datasets
	cout << "Cross Validation Index ... " << endl;
	vector <int> CVIndex;
	getCVIndex(cPar.prop, cPar.seed, sqcNAIndex, pheno, CVIndex);

	int NInitial = pheno.size();
	int NTraining = 0;
	vector<double> phenoTraining;
	for (int i = 0; i < NInitial; i++) {

		if (CVIndex[i] == 1) {

			phenoTraining.push_back(atof(pheno[i][0].c_str()));
			NTraining++;
		}
	}
	VectorXd phenotype = VectorXd::Map(&phenoTraining[0], NTraining);
	VectorXd phenotypeStd = standardizeVec(phenotype);

	//Summary data
	cout << "Summary data ... " << endl;
	string bgenEnd = "_v2.bgen";
	string bgen = cPar.bgenStart + cPar.chr + bgenEnd; 
	string summEnd = ".txt";
	string summOutput = cPar.outpath + "/" + cPar.outfile + "_chr" + cPar.chr + "_cv" + cPar.cv + summEnd;
	ofstream summFout(summOutput.c_str());
	MS.streamBgen(cPar.chr, NTraining, bgen, phenotypeStd, CVIndex,
						cPar.bgenMinMAF, cPar.bgenMinINFO, cPar.bgenMinHWE, cPar.bgenMinCallingRate, cPar.threads, 
						summFout);
	cout << "Summary data is OK!" << endl;
	return;
}
