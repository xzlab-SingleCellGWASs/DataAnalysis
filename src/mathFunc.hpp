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


#ifndef __MATHFUNC_H__
#define __MATHFUNC_H__

#include <eigen3/Eigen/Dense>
#include <vector>
#include <math.h>
#include <numeric>

using namespace std;
using namespace Eigen;

double CalcHWE(const int n_hom1, const int n_hom2, const int n_ab);

VectorXd standardizeVec(VectorXd x);

VectorXd standardizeVec(VectorXd x, double mean, double sd);

double dotProd(VectorXd x, VectorXd y);

bool SortCol2(const vector<double>& v1, const vector<double>& v2);

bool SortCol1(const vector<double>& v1, const vector<double>& v2);

double calcSd(VectorXd x);

double calcMSE(VectorXd y, VectorXd yHat);

double calcCOR(VectorXd x, VectorXd y);

#endif
