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


#ifndef __SNPEFF_H__
#define __SNPEFF_H__

#include <cstring>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <sys/time.h>
#include <math.h>
#include <stdlib.h>

#include "zlib.h"
#include "omp.h"
#include <boost/math/distributions/students_t.hpp>
#include <boost/math/distributions/chi_squared.hpp>

#include "Types.hpp"
#include "io.hpp"
#include "mathFunc.hpp"
#include "lmFunc.hpp"

using namespace std;
using namespace Eigen;

class MKSUMM {

public:

	string getMargin(int chrom, int physpos, string snpName, string allele1, string allele0,
						double reffrq, double info, uint64 filepos, int sampleSize, double dosage,
						vector <double> xVec, VectorXd y);
	
	string getSnpStatsBgen(uchar *buf, uint bufLen, const uchar *zBuf, uint zBufLen, uint Nbgen, int NTraining,
							const string &snpName, int chrom, int physpos, string allele1, string allele0, uint64 filepos,
							VectorXd phenotype, vector <int> CVIndex, double *lut,
							float bgenMinMAF, float bgenMinINFO, float bgenMinHWE, float bgenMinCall);

	void streamBgen(string chr, int NTraining, const string &bgen, VectorXd phenotype,
					vector <int> CVIndex, float bgenMinMAF, float bgenMinINFO, float bgenMinHWE, float bgenMinCall, int threads,
					ofstream &foutSummary);
};

class STEPSEL {

public:

	int snpSelChr(vector <vector <double> > snpInfo, float Len, vector <vector <double> > &snpInfoIdp);

	int getSnpStatsBgen(uchar *buf, uint bufLen, const uchar *zBuf, uint zBufLen, double *lut,
						int N, vector <int> CVIndex, string label, vector<double> &genotypeVec);

	int streamBgenChr(vector <vector <double> > summInfoChr, int chr, int N,
						vector <int> CVIndex, string begin, string end, int threads, string label,
						vector <vector <double> > &genoChr);

	int streamBgenSingle(vector <double> summInfoSingle, vector <int> CVIndex, string begin, string end, string label,
						vector <double> &genoVec);

	int getGenotype1(vector<vector <string> > summData, float Len, float threshold1, float threshold2,
					int NTraining, int NTesting, vector<int> CVIndex, string begin, string end, int threads,
					string label, MatrixXd &genoTotMat, vector<vector <double> > &summInfo);

	int getGenotype2(vector<vector <string> > summData, float Len, float threshold1, float threshold2,
					int NTraining, int NTesting, vector<int> CVIndex, string begin, string end, int threads,
					string label, VectorXd resid, vector <double> &pTot, vector<vector <double> > &snpIdp);
};

#endif
