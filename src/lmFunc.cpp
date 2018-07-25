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
#include <vector>
#include <string>
#include <iostream>
#include <set>
#include <algorithm>
#include <numeric>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <boost/math/distributions/students_t.hpp>
#include <boost/math/distributions/normal.hpp>

#include "mathFunc.hpp"
#include "snpEff.hpp"
#include "lmFunc.hpp"

using namespace std;
using namespace Eigen;

VectorXd mvlmResid(MatrixXd X, VectorXd y) 
{
	//* X and y is standardized.
	uint sampleSize = X.rows();
	VectorXd beta = X.transpose() * y / sampleSize;
	VectorXd epsilon = y - X * beta;

	return epsilon;
}

double lmTest(VectorXd X, VectorXd y)
{

	//* X and y is standardized.
	uint sampleSize = X.size();
	VectorXd yStd = standardizeVec(y);
	double z = X.dot(yStd) / sqrt(sampleSize);
	boost::math::students_t tDis(sampleSize);
	double p = 2 * cdf(complement(tDis, fabs(z)));
	
	return p;
}

int lmStepwise(MatrixXd XInit, VectorXd y, vector <vector <double> > summInfo, vector<int> CVIndex,
				double threshold, string begin, string end, vector <vector <double> > p, 
				vector <vector <double> > &summInfoSig, MatrixXd &XIter)
{

	STEPSEL SS;
	/*******Get snps (<threshold2) of multiple regression model*******/
	//*XInit and y is standardized.
	int sampleSize = XInit.rows();
	boost::math::normal_distribution<> nd(0.0, 1.0);
	XIter = XInit;
	
	vector <int> idx;
	for (int i = 0; i < summInfo.size(); i++) {

		// get the genotype of threshold2
		vector <double>  genoInd;
		SS.streamBgenSingle(summInfo[p[i][0]], CVIndex, begin, end, "Training", genoInd);
		VectorXd genoIndVec = standardizeVec(VectorXd::Map(&genoInd[0], genoInd.size()));

		VectorXd epsilonVec0 = mvlmResid(XIter, y);
		VectorXd epsilonVec1 = mvlmResid(XIter, genoIndVec);
		double pr = calcCOR(epsilonVec0, epsilonVec1);
		double z = 0.5 * log((1 + pr) / (1 - pr));
		double rho = 1 / sqrt(sampleSize - XIter.cols() - 2);
		double zeta = 0.5 * log((1 + rho) / (1 - rho));

		double pMSE = 1 - cdf(nd, sqrt(sampleSize - XIter.cols() - 3) * (abs(z) - zeta)) + cdf(nd, sqrt(sampleSize - XIter.cols() - 3) * (-abs(z) - zeta));
		cout << i+1 << "th snp, "<< "pMSE: " << pMSE << endl;
		if (pMSE > threshold) {

			break;
		}
		else {

			MatrixXd XMulti(sampleSize, XIter.cols() + 1);
			XMulti << XIter, genoIndVec;
			XIter = XMulti;
			idx.push_back((int) p[i][0]);
		}
	}
	sort(idx.begin(), idx.end());

	for (int i = 0; i < idx.size(); i++) {

		summInfoSig.push_back(summInfo[idx[i]]);
	}
	return 0;
}
