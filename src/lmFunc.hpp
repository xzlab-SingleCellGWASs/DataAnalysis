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

#ifndef __LMFUNC_H__
#define __LMFUNC_H__

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

#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/normal.hpp>

#include "mathFunc.hpp"
#include "snpEff.hpp"

using namespace std;
using namespace Eigen;

VectorXd mvlmResid(MatrixXd X, VectorXd y);

double lmTest(VectorXd X, VectorXd y);

int lmStepwise(MatrixXd XInit, VectorXd y, vector <vector <double> > SummInfo, vector<int> CVIndex,
	double threshold, string begin, string end, vector <vector <double> > p,
	vector <vector <double> > &SummInfoSig, MatrixXd &XIter);

#endif
