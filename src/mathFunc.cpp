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
#include <math.h>
#include <numeric>
#include <iostream>

#include "mathFunc.hpp"

using namespace std;
using namespace Eigen;

double CalcHWE(const int n_hom1, const int n_hom2, const int n_ab) {

	// "AA" is the rare allele.
	int n_aa = n_hom1 < n_hom2 ? n_hom1 : n_hom2;
	int n_bb = n_hom1 < n_hom2 ? n_hom2 : n_hom1;

	int rare_copies = 2 * n_aa + n_ab;
	int genotypes = n_ab + n_bb + n_aa;

	double *het_probs = (double *)malloc((rare_copies + 1) * sizeof(double));
	//if (het_probs == NULL)
	//	cerr << "Internal error: SNP-HWE: Unable to allocate array" << endl;

	int i;
	for (i = 0; i <= rare_copies; i++)
		het_probs[i] = 0.0;

	// Start at midpoint.
	// XZ modified to add (long int)
	int mid = ((long int)rare_copies *
		(2 * (long int)genotypes - (long int)rare_copies)) /
		(2 * (long int)genotypes);

	// Check to ensure that midpoint and rare alleles have same
	// parity.
	if ((rare_copies & 1) ^ (mid & 1))
		mid++;
	
	int curr_hets = mid;
	int curr_homr = (rare_copies - mid) / 2;
	int curr_homc = genotypes - curr_hets - curr_homr;

	het_probs[mid] = 1.0;
	double sum = het_probs[mid];
	for (curr_hets = mid; curr_hets > 1; curr_hets -= 2) {
		het_probs[curr_hets - 2] = het_probs[curr_hets] * curr_hets *
			(curr_hets - 1.0) /
			(4.0 * (curr_homr + 1.0) * (curr_homc + 1.0));
		sum += het_probs[curr_hets - 2];

		// Two fewer heterozygotes for next iteration; add one
		// rare, one common homozygote.
		curr_homr++;
		curr_homc++;
	}

	curr_hets = mid;
	curr_homr = (rare_copies - mid) / 2;
	curr_homc = genotypes - curr_hets - curr_homr;
	for (curr_hets = mid; curr_hets <= rare_copies - 2; curr_hets += 2) {
		het_probs[curr_hets + 2] = het_probs[curr_hets] * 4.0 * curr_homr *
			curr_homc /
			((curr_hets + 2.0) * (curr_hets + 1.0));
		sum += het_probs[curr_hets + 2];

		// Add 2 heterozygotes for next iteration; subtract
		// one rare, one common homozygote.
		curr_homr--;
		curr_homc--;
	}

	for (i = 0; i <= rare_copies; i++)
		het_probs[i] /= sum;

	double p_hwe = 0.0;

	// p-value calculation for p_hwe.
	for (i = 0; i <= rare_copies; i++) {
		if (het_probs[i] > het_probs[n_ab])
			continue;
		p_hwe += het_probs[i];
	}

	p_hwe = p_hwe > 1.0 ? 1.0 : p_hwe;


	free(het_probs);

	return p_hwe;
}


VectorXd standardizeVec(VectorXd x) {

	int n = x.size();
	VectorXd xMean(x);	
	xMean.fill(x.mean());
	double xMulti = (x - xMean).dot(x - xMean);
	double xSd = sqrt(xMulti / n);
	VectorXd xStandard = (x - xMean) / xSd;

	return xStandard;
}

VectorXd standardizeVec(VectorXd x, double mean, double sd) {

	int n = x.size();
	VectorXd xMean(n);
	xMean.fill(mean);
	VectorXd xStandard = (x - xMean) / sd;

	return xStandard;
}


double dotProd(VectorXd x, VectorXd y)
{
	
	int N = x.size();
	double dot = x.transpose() * y;
	double dotN = dot * N;
	return dotN;
}

bool SortCol2(const vector<double>& v1, const vector<double>& v2)
{

	return v1[2] < v2[2];
}

bool SortCol1(const vector<double>& v1, const vector<double>& v2)
{

	return v1[1] < v2[1];
}

double calcSd(VectorXd x)
{
	
	double ans = 0;
	double mean = x.mean();
	for (int i = 0; i < x.size(); i++) {

		ans += (x(i) - mean) * (x(i) - mean);
	}
	double sd = sqrt(ans / (x.size()-1));
	
	return sd;
}

double calcMSE(VectorXd y, VectorXd yHat)
{

	double ss = dotProd(y - yHat, y - yHat);
	double MSE = ss / pow(y.size(), 2);
	return MSE;
}

double calcCOR(VectorXd x, VectorXd y)
{

	int len = x.size();
	VectorXd xMean(len), yMean(len);
	xMean.fill(x.mean());
	yMean.fill(y.mean());
	double cov = (y - yMean).dot(x - xMean);
	double xVar = (x - xMean).dot(x - xMean);
	double yVar = (y - yMean).dot(y - yMean);
	double COR = cov / sqrt(xVar * yVar);
	return COR;
}