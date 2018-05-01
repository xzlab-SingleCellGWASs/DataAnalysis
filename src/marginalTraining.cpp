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

#include "Types.hpp"
#include "io.hpp"
#include "marginalTraining.hpp"
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>


using namespace std;
using namespace Eigen;

int getMarginTraining(int chrom, int physpos, string allele1, string allele0, double maf, double info, double filepos,
	VectorXd x, VectorXd y, MatrixXd sqc, vector <double> &snpMarginalInfo)
{

	int sampleSize = x.rows();

	/***Parameter Estimation***/
	MatrixXd X(x.rows(), 1 + sqc.cols());
	X << x, sqc;
	MatrixXd XXInv = (X.transpose() * X).inverse();
	VectorXd beta = XXInv * X.transpose() * y;

	/***Hypothesis Testing (t test)***/
	VectorXd sigma = y - X * beta;
	double sigmaSum = (sigma.transpose()*sigma).sum();
	double sigmaX = XXInv.diagonal()(0) * sigmaSum;
	int df = sampleSize - X.cols();
	double se = sqrt(sigmaX / df);
	double t = beta(0) / se;
	boost::math::students_t dist(df);
	double p = 2 * cdf(complement(dist, fabs(t)));

	/***Output***/
	snpMarginalInfo.resize(9);
	snpMarginalInfo[0] = chrom;
	snpMarginalInfo[1] = physpos;
	snpMarginalInfo[2] = baseToCode(allele1);
	snpMarginalInfo[3] = baseToCode(allele0);
	snpMarginalInfo[4] = maf;
	snpMarginalInfo[5] = info;
	snpMarginalInfo[6] = p;
	snpMarginalInfo[7] = t;
	snpMarginalInfo[8] = (double)filepos;

	return 0;
}

int getSnpStatsBgenTraining(uchar *buf, uint bufLen, const uchar *zBuf, uint zBufLen, uint Nbgen, int NTraining,
	const string &snpName, int chrom, int physpos, string allele1, string allele0, int64 filepos,
	VectorXd phenotype, MatrixXd sqc, vector <int> CVIndex,
	float bgenMinMAF, float bgenMinINFO, float bgenMinHW, float bgenCallingRate,
	vector <double> &snpMarginalInfo)
{

	/********** decompress and check genotype probability block **********/
	uLongf destLen = bufLen;
	if (uncompress(buf, &destLen, zBuf, zBufLen) != Z_OK || destLen != bufLen) {
		cerr << "ERROR: uncompress() failed" << endl;
		exit(1);
	}
	uchar *bufAt = buf;
	uint N = bufAt[0] | (bufAt[1] << 8) | (bufAt[2] << 16) | (bufAt[3] << 24); bufAt += 4;
	if (N != Nbgen) {
		cerr << "ERROR: " << snpName << " has N = " << N << " (mismatch with header block)" << endl;
		exit(1);
	}
	uint K = bufAt[0] | (bufAt[1] << 8); bufAt += 2;
	if (K != 2U) {
		cerr << "ERROR: " << snpName << " has K = " << K << " (non-bi-allelic)" << endl;
		exit(1);
	}
	uint Pmin = *bufAt; bufAt++;
	if (Pmin != 2U) {
		cerr << "ERROR: " << snpName << " has minimum ploidy = " << Pmin << " (not 2)" << endl;
		exit(1);
	}
	uint Pmax = *bufAt; bufAt++;
	if (Pmax != 2U) {
		cerr << "ERROR: " << snpName << " has maximum ploidy = " << Pmax << " (not 2)" << endl;
		exit(1);
	}
	for (uint i = 0; i < N; i++) {
		uint ploidyMiss = *bufAt; bufAt++;
		if (ploidyMiss != 2U) {
			cerr << "ERROR: " << snpName << " has ploidy/missingness byte = " << ploidyMiss
				<< " (not 2)" << endl;
			exit(1);
		}
	}
	uint Phased = *bufAt; bufAt++;
	if (Phased != 0U) {
		cerr << "ERROR: " << snpName << " has Phased = " << Phased << " (not 0)" << endl;
		exit(1);
	}
	uint B = *bufAt; bufAt++;
	if (B != 8U) {
		cerr << "ERROR: " << snpName << " has B = " << B << " (not 8)" << endl;
		exit(1);
	}

	/********** compute and filter MAF and INFO **********/

	snpMarginalInfo.resize(9);
	//MAF and INFO
	double lut[256];
	for (int i = 0; i <= 255; i++)
		lut[i] = i / 255.0;

	double sum_eij = 0, sum_fij_minus_eij2 = 0;
	vector<double> genotype; genotype.reserve(N);
	vector<double> P11; P11.reserve(N);
	vector<double> P10; P10.reserve(N);
	vector<double> P00;	P00.reserve(N);
	for (uint i = 0; i < N; i++) {

		double p11 = lut[*bufAt]; bufAt++;
		double p10 = lut[*bufAt]; bufAt++;
		double p00 = 1 - p11 - p10;
		double eij = 2 * p11 + p10;
		genotype.push_back(2 * p00 + p10);
		P11.push_back(p11);
		P10.push_back(p10);
		P00.push_back(p00);
		double fij = 4 * p11 + p10;
		sum_eij += eij;
		sum_fij_minus_eij2 += fij - eij * eij;
	}

	double reffrq = sum_eij / (2 * N);
	double info = reffrq == 0 || reffrq == 1 ? 1 :
		1 - sum_fij_minus_eij2 / (2 * N * reffrq * (1 - reffrq));
	double maf = min(reffrq, 1 - reffrq);

	//Selection by MAF and INFO
	if (maf > bgenMinMAF && info > bgenMinINFO) {

		//HWE
		double P11_sum = accumulate(P11.begin(), P11.end(), 0.0); //cout << "snpName: " << snpName << "P11_sum: " << P11_sum << endl;
		double P10_sum = accumulate(P10.begin(), P10.end(), 0.0); //cout << "snpName: " << snpName << "P10_sum: " << P10_sum << endl;
		double P00_sum = accumulate(P00.begin(), P00.end(), 0.0); //cout << "snpName: " << snpName << "P00_sum: " << P00_sum << endl;

		double p = (2 * P11_sum + P10_sum) / (2 * (P11_sum + P10_sum + P00_sum));
		double q = 1 - p;
		double e11 = pow(p, 2) * N;
		double e10 = 2 * p * q * N;
		double e00 = pow(q, 2) * N;

		double chi_sq = pow((P11_sum - e11), 2) / e11 + pow((P10_sum - e10), 2) / e10 + pow((P00_sum - e00), 2) / e00; //cout << "chisq: " << chi_sq << endl;
		boost::math::chi_squared dist(1);
		double chi_p = 1 - cdf(dist, chi_sq);

		//Selection by HWE
		if (chi_p >= bgenMinHW) {

			vector<double> genoTraining;
			vector<double> P11Training;
			vector<double> P10Training;
			vector<double> P00Training;
			double dosageSum = 0;
			for (int i = 0; i < N; i++) {

				if (CVIndex[i] == 1) {

					genoTraining.push_back(genotype[i]);
					P11Training.push_back(P11[i]);
					P10Training.push_back(P10[i]);
					P00Training.push_back(P00[i]);
					dosageSum += genotype[i];
				}
			}
			int genoTrainingSize = genoTraining.size();
			int noncall = 0;
			for (int i = 0; i < genoTrainingSize; i++) {

				if (P11Training[i] == 0 && P10Training[i] == 0 && P00Training[i] == 0)  noncall++;
			}
			double callingRate = 1 - (noncall * 0.1) / (genoTrainingSize * 0.1);

			//Selection by calling rate
			if (callingRate > bgenCallingRate) {

				VectorXd genoTrainingVec = VectorXd::Map(&genoTraining[0], genoTraining.size());
				getMarginTraining(chrom, physpos, allele1, allele0, reffrq, info, filepos,
					genoTrainingVec, phenotype, sqc, snpMarginalInfo);

			}
			else {

				fill(snpMarginalInfo.begin(), snpMarginalInfo.end(), 99);
			}
		}
		else {

			fill(snpMarginalInfo.begin(), snpMarginalInfo.end(), 99);
		}
	}
	else {

		fill(snpMarginalInfo.begin(), snpMarginalInfo.end(), 99);
	}
	return 0;
}

bool SortCol2(const vector<double>& v1, const vector<double>& v2)
{

	return v1[2] < v2[2];
}

int snpSel(vector <vector <double> > snpInfo, int Len,
	vector <vector <double> > &snpInfoIdp)
{

	sort(snpInfo.begin(), snpInfo.end(), SortCol2);

	int snpNum = snpInfo.size();
	vector <int> b(snpNum);
	fill(b.begin(), b.end(), 0);

	for (int i = 0; i < snpNum; i++) {

		for (int j = i + 1; j < snpNum; j++) {

			int dis = abs(snpInfo[i][1] - snpInfo[j][1]);
			if (dis < Len) {

				b[j] += 1;
			}
		}
	}

	for (int i = 0; i < snpNum; i++) {

		if (b[i] == 0) snpInfoIdp.push_back(snpInfo[i]);
	}

	return 0;
}

int streamBgenTraining(string chr, int NTraining, const string &bgen, VectorXd phenotype,
	MatrixXd sqc, vector <int> CVIndex, double bgenMinMAF, double bgenMinINFO,
	double bgenMinHW, double bgenCallingRate, int threads,
	float threshold1, float threshold2,
	ofstream &foutSummary, ofstream &foutInfoI, ofstream &foutInfoII)
{

	//Times
	struct timeval start, end;
	//Input
	FILE* fin_x = fopen(bgen.c_str(), "rb");

	/********** READ HEADER **********/
	uint offset; fread(&offset, 4, 1, fin_x); //cout << "offset: " << offset << endl;
	uint L_H; fread(&L_H, 4, 1, fin_x); //cout << "L_H: " << L_H << endl;
	uint Mbgen; fread(&Mbgen, 4, 1, fin_x);  cout << "snpBlocks (Mbgen): " << Mbgen << endl;
	uint Nbgen; fread(&Nbgen, 4, 1, fin_x);  cout << "samples (Nbgen): " << Nbgen << endl;

	char magic[5]; fread(magic, 1, 4, fin_x); magic[4] = '\0'; //cout << "magic bytes: " << string(magic) << endl;
	fseek(fin_x, L_H - 20, SEEK_CUR); //cout << "skipping L_H-20 = " << L_H-20 << " bytes (free data area)" << endl;
	uint flags; fread(&flags, 4, 1, fin_x); //cout << "flags: " << flags << endl;
	uint CompressedSNPBlocks = flags & 3;  //cout << "CompressedSNPBlocks: " << CompressedSNPBlocks << endl;
	assert(CompressedSNPBlocks == 1); // REQUIRE CompressedSNPBlocks==1
	uint Layout = (flags >> 2) & 0xf;  //cout << "Layout: " << Layout << endl;
	assert(Layout == 1 || Layout == 2); // REQUIRE Layout==1 or Layout==2

	fseek(fin_x, offset + 4, SEEK_SET);
	int fInitial = ftell(fin_x);

	/********** READ SNP BLOCKS IN BATCHES **********/

	const int B_MAX = 400; // number of SNPs to process in one batch (for multi-threading)

	char snpID[65536], rsID[65536], chrStr[65536];
	char *allele1, *allele0;
	uint maxLA = 65536, maxLB = 65536;
	allele1 = (char *)malloc(maxLA + 1);
	allele0 = (char *)malloc(maxLB + 1);

	// during single-threaded reading of block, store SNP data for later multi-threaded processing
	vector < vector <uchar> > bufs(threads);
	vector <string> snpNames(B_MAX);
	vector <int> chroms(B_MAX);
	vector <int> bps(B_MAX);
	vector <int64> fps(B_MAX);
	vector <string> allele1s(B_MAX), allele0s(B_MAX);
	vector <vector <double> > snpMarginalInfoBlock(B_MAX);
	vector <vector <double> > snpMarginalInfo;
	vector <vector <double> > snpMarginalInfoI;
	vector <vector <double> > snpMarginalInfoII;
	vector < vector <uchar> > zBufs(B_MAX);
	vector <uint> zBufLens(B_MAX), bufLens(B_MAX);

	int B = 0;

	gettimeofday(&start, NULL);
	//each snp  Mbgen
	for (uint mbgen = 0; mbgen < Mbgen; mbgen++) {

		fps[B] = ftell(fin_x) - fInitial;

		ushort LS; fread(&LS, 2, 1, fin_x);  //cout << "LS: " << LS << " " << std::flush;
		fread(snpID, 1, LS, fin_x); snpID[LS] = '\0'; // cout << "snpID: " << string(snpID) << " " << std::flush;
		ushort LR; fread(&LR, 2, 1, fin_x); // cout << "LR: " << LR << " " << std::flush;
		fread(rsID, 1, LR, fin_x); rsID[LR] = '\0';  //cout << "rsID: " << string(rsID) << " " << std::flush;
		snpNames[B] = string(rsID) == "." ? snpID : rsID;

		ushort LC; fread(&LC, 2, 1, fin_x); // cout << "LC: " << LC << " " << std::flush;
		fread(chrStr, 1, LC, fin_x); chrStr[LC] = '\0';
		int chrom;
		if (sscanf(chrStr, "%d", &chrom) != 1 || !(1 <= chrom && chrom <= 22)) {

			cerr << "ERROR: Invalid chrom (expecting integer 1-22): " << string(chrStr) << endl;
			exit(1);
		}
		chroms[B] = chrom;

		uint physpos; fread(&physpos, 4, 1, fin_x); // cout << "physpos: " << physpos << " " << std::flush;
		bps[B] = physpos;

		ushort K; fread(&K, 2, 1, fin_x); //cout << "K: " << K << endl;
		if (K != 2) {

			cerr << "ERROR: Non-bi-allelic variant found: " << K << " alleles" << endl;
			exit(1);
		}

		uint LA; fread(&LA, 4, 1, fin_x);  //cout << "LA: " << LA << " " << std::flush;
		if (LA > maxLA) {

			maxLA = 2 * LA;
			free(allele1);
			allele1 = (char *)malloc(maxLA + 1);
		}
		fread(allele1, 1, LA, fin_x); allele1[LA] = '\0';
		allele1s[B] = string(allele1);

		uint LB; fread(&LB, 4, 1, fin_x);  //cout << "LB: " << LB << " " << std::flush;
		if (LB > maxLB) {

			maxLB = 2 * LB;
			free(allele0);
			allele0 = (char *)malloc(maxLB + 1);
		}
		fread(allele0, 1, LB, fin_x); allele0[LB] = '\0';
		allele0s[B] = string(allele0);

		uint C; fread(&C, 4, 1, fin_x); //cout << "C: " << C << endl;
		if (C > zBufs[B].size()) zBufs[B].resize(C - 4);
		uint D; fread(&D, 4, 1, fin_x); //cout << "D: " << D << endl;
		zBufLens[B] = C - 4; bufLens[B] = D;

		fread(&zBufs[B][0], 1, C - 4, fin_x);

		B++;

		if (B == B_MAX || mbgen + 1 == Mbgen) { // process the block of SNPs using multi-threading

			omp_set_num_threads(threads);

#pragma omp parallel for 
			for (int b = 0; b < B; b++) {

				int t = omp_get_thread_num();
				if (bufLens[b] > bufs[t].size()) bufs[t].resize(bufLens[b]);
				vector <double> snpMarginalInfo;
				getSnpStatsBgenTraining(&bufs[t][0], bufLens[b], &zBufs[b][0], zBufLens[b], Nbgen, NTraining,
					snpNames[b], chroms[b], bps[b], allele1s[b], allele0s[b], fps[b],
					phenotype, sqc, CVIndex,
					bgenMinMAF, bgenMinINFO, bgenMinHW, bgenCallingRate, snpMarginalInfo);
				int snpMarginalInfoSize = snpMarginalInfo.size();
				snpMarginalInfoBlock[b].resize(snpMarginalInfoSize);
				snpMarginalInfoBlock[b] = snpMarginalInfo;
			}

			/***Output the summary data (MAF, INFO and HWE)***/
			for (int b = 0; b < B; b++) {

				if (snpMarginalInfoBlock[b][0] != 99) {

					string standardOutput = standardOutputSummaryData(snpMarginalInfoBlock[b], snpNames[b]);
					foutSummary << standardOutput;
					snpMarginalInfo.push_back(snpMarginalInfoBlock[b]);
				}
			}

			B = 0; // reset current block size
		}

		if (mbgen % 250000 == 249999) {

			gettimeofday(&end, NULL);
			int timeuse = (end.tv_sec + 1e-6 * end.tv_usec) - (start.tv_sec + 1e-6 * start.tv_usec);
			cout << "mbgen: " << mbgen + 1 << "\t" << "Block time: " << timeuse << endl;
		}
	}

	cout << "After flitering by MAF, INFO and HEW, chrom " << chr << " includes " << snpMarginalInfo.size() << " SNPs." << endl;

	for (int i = 0; i < snpMarginalInfo.size(); i++) {

		if (snpMarginalInfo[i][6] < threshold1 && snpMarginalInfo[i][6] > threshold2) {

			vector <double> snpMarginalFliter(4);
			snpMarginalFliter[0] = snpMarginalInfo[i][0];
			snpMarginalFliter[1] = snpMarginalInfo[i][1];
			snpMarginalFliter[2] = snpMarginalInfo[i][6];
			snpMarginalFliter[3] = snpMarginalInfo[i][8];
			snpMarginalInfoI.push_back(snpMarginalFliter);
		}

		if (snpMarginalInfo[i][6] < threshold2) {

			vector <double> snpMarginalFliter(4);
			snpMarginalFliter[0] = snpMarginalInfo[i][0];
			snpMarginalFliter[1] = snpMarginalInfo[i][1];
			snpMarginalFliter[2] = snpMarginalInfo[i][6];
			snpMarginalFliter[3] = snpMarginalInfo[i][8];
			snpMarginalInfoII.push_back(snpMarginalFliter);
		}

	}

	cout << "Chrom " << chr << " includes " << snpMarginalInfoI.size() << " SNPs > " << threshold2 << " and < " << threshold1 << "." << endl;
	cout << "Chrom " << chr << " includes " << snpMarginalInfoII.size() << " SNPs < " << threshold2 << "." << endl;

	vector <vector <double> > snpMarginalInfoIdpI;
	vector <vector <double> > snpMarginalInfoIdpII;
	snpSel(snpMarginalInfoI, 1e6, snpMarginalInfoIdpI);
	snpSel(snpMarginalInfoII, 1e6, snpMarginalInfoIdpII);

	for (int i = 0; i < snpMarginalInfoIdpI.size(); i++) foutInfoI << standardOutputInfo(snpMarginalInfoIdpI[i]);
	for (int i = 0; i < snpMarginalInfoIdpII.size(); i++) foutInfoII << standardOutputInfo(snpMarginalInfoIdpII[i]);

	cout << "Chrom " << chr << " includes " << snpMarginalInfoIdpI.size() << " indpendent SNPs > " << threshold2 << " and < " << threshold1 << "." << endl;
	cout << "Chrom " << chr << " includes " << snpMarginalInfoIdpII.size() << " indpendent SNPs < " << threshold2 << "." << endl;

	free(allele0);
	free(allele1);
	return 0;
}