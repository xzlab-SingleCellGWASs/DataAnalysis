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


#ifndef __MARGINALTRAINING_H__
#define __MARGINALTRAINING_H__

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


using namespace std;
using namespace Eigen;

int getMarginTraining(int chrom, int physpos, string allele1, string allele0, double maf, double info, double filepos,
						VectorXd x, VectorXd y, MatrixXd sqc, 
						vector <double> &snpMarginalInfo);

int getSnpStatsBgenTraining(uchar *buf, uint bufLen, const uchar *zBuf, uint zBufLen, uint Nbgen, int NTraining,
							const string &snpName, int chrom, int physpos, string allele1, string allele0, double filepos,
							VectorXd phenotype, MatrixXd sqc, vector <int> CVIndex,
							float bgenMinMAF, float bgenMinINFO, float bgenMinHW, float bgenCallingRate,
							vector <double> &snpMarginalInfo);

bool SortCol2(const vector<double>& v1, const vector<double>& v2);

int snpSel(vector <vector <double> > snpInfo, int Len, vector <vector <double> > &snpInfoIdp);

int streamBgenTraining(string chr, int NTraining, const string &bgen, VectorXd phenotype,
						MatrixXd sqc, vector <int> CVIndex, double bgenMinMAF, double bgenMinINFO,
						double bgenMinHW, double bgenCallingRate, int threads,
						float threshold1, float threshold2,
						ofstream &foutSummary, ofstream &foutInfoI, ofstream &foutInfoII);

#endif
