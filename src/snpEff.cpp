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
#include "snpEff.hpp"

using namespace std;
using namespace Eigen;

string MKSUMM::getMargin(int chrom, int physpos, string snpName, string allele1, string allele0,
	double reffrq, double info, uint64 filepos, int sampleSize, double dosage,
	vector <double> xVec, VectorXd y)
{

	/***Transformation***/
	VectorXd x = VectorXd::Map(&xVec[0], sampleSize);
	VectorXd xStd = standardizeVec(x);
	double xMean = x.mean();
	double xSd = calcSd(x);
	double xty = xStd.transpose() * y;
	double beta = xty / sampleSize;
	boost::math::students_t tDis(sampleSize);
	double t = sqrt(sampleSize) * beta;
	try{

		double p = 2 * cdf(complement(tDis, fabs(t)));
		string chr = "chr" + to_string(chrom);

		/***Output***/
		ostringstream ossSumm;
		ossSumm << chr << "\t" << physpos << "\t" << allele1 << "\t" << allele0 << "\t" << reffrq <<
			"\t" << info << "\t" << snpName << "\t" << dosage << "\t" << p << "\t" << beta << "\t" << filepos << 
			"\t" << exp(beta) << "\t" << xMean << "\t" << xSd;
		ossSumm << endl;

		return ossSumm.str();
	}
	catch (exception& e) {
		return "";
	}
}

string MKSUMM::getSnpStatsBgen(uchar *buf, uint bufLen, const uchar *zBuf, uint zBufLen, uint Nbgen, int NTraining,
	const string &snpName, int chrom, int physpos, string allele1, string allele0, uint64 filepos,
	VectorXd phenotype, vector <int> CVIndex, double *lut,
	float bgenMinMAF, float bgenMinINFO, float bgenMinHWE, float bgenMinCall)
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

	/********** compute and filter MAF, INFO and calling rate **********/
	//MAF and INFO
	double sum_eij = 0, sum_fij_minus_eij2 = 0;
	vector<double> genotype; genotype.reserve(NTraining);
	double p11_sum = 0, p10_sum = 0, p00_sum = 0;
	double dosage = 0;
	double nCall = 0;
	for (uint i = 0; i < N; i++) {

		double p11 = lut[*bufAt]; bufAt++;
		double p10 = lut[*bufAt]; bufAt++;
		double p00 = 1 - p11 - p10;
		double eij = 2 * p11 + p10;
		double fij = 4 * p11 + p10;
		sum_eij += eij;
		sum_fij_minus_eij2 += fij - eij * eij;

		if (CVIndex[i] == 1) {

			double aij = 2 * p00 + p10;
			dosage += aij;
			p11_sum += p11;
			p10_sum += p10;
			p00_sum += p00;
			genotype.push_back(aij);
			if (p11 == 0 || p10 == 0 || p00 == 0) nCall++;
		}
	}

	double reffrq = sum_eij / (2 * N);
	double info = reffrq == 0 || reffrq == 1 ? 1 : 1 - sum_fij_minus_eij2 / (2 * N * reffrq * (1 - reffrq));
	double maf = min(reffrq, 1 - reffrq);
	double callRate = (nCall * 0.1) / (NTraining * 0.1);

	if (maf < bgenMinMAF || info < bgenMinINFO || callRate < bgenMinCall) return "";
	
	double hwe_p = CalcHWE(ceil(p11_sum), ceil(p00_sum), ceil(p10_sum));

	if (hwe_p < bgenMinHWE) return "";
	//cout << "snpName: " << snpName << " maf: " << maf << " hwe_p: " << hwe_p << endl;
	//training genotype
	ostringstream ossSumm;
	ossSumm << getMargin(chrom, physpos, snpName, allele1, allele0, reffrq, info, filepos, NTraining, dosage, genotype, phenotype);

	return ossSumm.str();
}

void MKSUMM::streamBgen(string chr, int NTraining, const string &bgen, VectorXd phenotype,
	vector <int> CVIndex, float bgenMinMAF, float bgenMinINFO, float bgenMinHWE, float bgenMinCall, int threads,
	ofstream &foutSummary)
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

	const int B_MAX = 1000; // number of SNPs to process in one batch (for multi-threading)

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
	vector <uint64> fps(B_MAX);
	vector <string> allele1s(B_MAX), allele0s(B_MAX);
	vector <string > snpMarginalInfoBlock(B_MAX);
	vector < vector <uchar> > zBufs(B_MAX);
	vector <uint> zBufLens(B_MAX), bufLens(B_MAX);

	int B = 0;

	double lut[256];
	for (int i = 0; i <= 255; i++)
		lut[i] = i / 255.0;

	gettimeofday(&start, NULL);
	//each snp  Mbgen
	for (uint mbgen = 0; mbgen < Mbgen; mbgen++) {
	
		fps[B] = ftell(fin_x) - fInitial; //cout << "fps[B]: " << fps[B] << endl;
		ushort LS; fread(&LS, 2, 1, fin_x);  //cout << "LS: " << LS << " " << std::flush;
		fread(snpID, 1, LS, fin_x); snpID[LS] = '\0'; // cout << "snpID: " << string(snpID) << " " << std::flush;
		ushort LR; fread(&LR, 2, 1, fin_x); // cout << "LR: " << LR << " " << std::flush;
		fread(rsID, 1, LR, fin_x); rsID[LR] = '\0';  //cout << "rsID: " << string(rsID) << " " << std::flush;
		snpNames[B] = string(rsID) == "." ? snpID : rsID; //cout << "snpNames[B]: " << snpNames[B] << "fps[B]: " << fps[B] << endl;

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
												//cout << "B: " << B <<  endl;
												//cout << mbgen << " ";
			omp_set_num_threads(threads);
#pragma omp parallel for schedule(dynamic)
			for (int b = 0; b < B; b++) {

				int t = omp_get_thread_num();
				if (bufLens[b] > bufs[t].size()) bufs[t].resize(bufLens[b]);

				snpMarginalInfoBlock[b] = getSnpStatsBgen(&bufs[t][0], bufLens[b], &zBufs[b][0], zBufLens[b],
						Nbgen, NTraining, snpNames[b], chroms[b], bps[b], allele1s[b], allele0s[b], fps[b],
						phenotype, CVIndex, lut,
						bgenMinMAF, bgenMinINFO, bgenMinHWE, bgenMinCall);
				//cout << chroms[b] << endl;
			}

			for (int b = 0; b < B; b++) foutSummary << snpMarginalInfoBlock[b];

			B = 0; // reset current block size
		}

		if (mbgen % 100000 == 99999) {

			gettimeofday(&end, NULL);
			int timeuse = (end.tv_sec + 1e-6 * end.tv_usec) - (start.tv_sec + 1e-6 * start.tv_usec);
			cout << "mbgen: " << mbgen + 1 << "\t" << "Block time: " << timeuse << endl;
		}
	}

	free(allele0);
	free(allele1);
	foutSummary.close();
	return;
}


int STEPSEL::snpSelChr(vector <vector <double> > snpInfo, float Len,
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

	sort(snpInfoIdp.begin(), snpInfoIdp.end(), SortCol1);
	return 0;
}

int STEPSEL::getSnpStatsBgen(uchar *buf, uint bufLen, const uchar *zBuf, uint zBufLen, double *lut,
	int N, vector <int> CVIndex, string label,	vector<double> &genotypeVec)
{

	/********** decompress and check genotype probability block **********/
	uLongf destLen = bufLen;
	uncompress(buf, &destLen, zBuf, zBufLen);
	uchar *bufAt = buf;
	bufAt += 10 + N;
	/********** get the genotype by the label **********/
	if (label == "Training") {

		for (uint i = 0; i < N; i++) {

			double p11 = lut[*bufAt]; bufAt++;
			double p10 = lut[*bufAt]; bufAt++;
			double p00 = 1 - p11 - p10;
			if (CVIndex[i] == 1) {

				genotypeVec.push_back(2 * p00 + p10);
			}
		}
	}
	if (label == "Testing") {

		for (uint i = 0; i < N; i++) {

			double p11 = lut[*bufAt]; bufAt++;
			double p10 = lut[*bufAt]; bufAt++;
			double p00 = 1 - p11 - p10;
			if (CVIndex[i] == 0) {

				genotypeVec.push_back(2 * p00 + p10);
			}
		}
	}
	return 0;
}

int STEPSEL::streamBgenChr(vector <vector <double> > summInfoChr, int chr, int N,
	vector <int> CVIndex, string begin, string end, int threads, string label,
	vector <vector <double> > &genoChr)
{

	string bgen = begin + to_string(chr) + end;
	FILE* fin_x = fopen(bgen.c_str(), "rb");
	/********** READ HEADER **********/
	uint offset; fread(&offset, 4, 1, fin_x);
	uint L_H; fread(&L_H, 4, 1, fin_x);
	uint Mbgen; fread(&Mbgen, 4, 1, fin_x);
	uint Nbgen; fread(&Nbgen, 4, 1, fin_x);

	char magic[5]; fread(magic, 1, 4, fin_x); magic[4] = '\0';
	fseek(fin_x, L_H - 20, SEEK_CUR);
	uint flags; fread(&flags, 4, 1, fin_x);

	/********** READ SNP BLOCKS IN BATCHES **********/
	const int B_MAX = 5; // number of SNPs to process in one batch (for multi-threading)
	char snpID[65536], rsID[65536], chrStr[65536];
	char *allele1, *allele0;
	uint maxLA = 65536, maxLB = 65536;
	allele1 = (char *)malloc(maxLA + 1);
	allele0 = (char *)malloc(maxLB + 1);

	double lut[256];
	for (int i = 0; i <= 255; i++)
		lut[i] = i / 255.0;

	// during single-threaded reading of block, store SNP data for later multi-threaded processing
	vector < vector <uchar> > bufs(threads);
	vector < vector <uchar> > zBufs(B_MAX);
	vector <uint> zBufLens(B_MAX), bufLens(B_MAX);

	int B = 0;

	//each snp  Mbgen
	int snpNum = summInfoChr.size();
	// cout << "snpNum: " << snpNum << endl;
	for (uint mbgen = 0; mbgen < snpNum; mbgen++) {

		fseek(fin_x, offset + 4, SEEK_SET);//cout << "start: " << ftell(fin_x) << endl;
		fseek(fin_x, summInfoChr[mbgen][3], SEEK_CUR);//cout << "current: " << ftell(fin_x) << endl;
		ushort LS; fread(&LS, 2, 1, fin_x);
		fread(snpID, 1, LS, fin_x);
		ushort LR; fread(&LR, 2, 1, fin_x);
		fread(rsID, 1, LR, fin_x);
		ushort LC; fread(&LC, 2, 1, fin_x);
		fread(chrStr, 1, LC, fin_x); chrStr[LC] = '\0';

		uint physpos; fread(&physpos, 4, 1, fin_x); //cout << "physpos: " << physpos << endl;
		ushort K; fread(&K, 2, 1, fin_x); //cout << "K: " << K << endl;
		uint LA; fread(&LA, 4, 1, fin_x); //cout << "LA: " << LA << endl;
		fread(allele1, 1, LA, fin_x); //cout << "allele1: " << allele1 << endl;
		uint LB; fread(&LB, 4, 1, fin_x); //cout << "LB: " << LB << endl;
		fread(allele0, 1, LB, fin_x); //cout << "allele0: " << allele0 << endl;

		uint C; fread(&C, 4, 1, fin_x); //cout << "C: " << C << endl;
		if (C > zBufs[B].size()) zBufs[B].resize(C - 4);
		uint D; fread(&D, 4, 1, fin_x); //cout << "D: " << D << endl;
		zBufLens[B] = C - 4; bufLens[B] = D;
		fread(&zBufs[B][0], 1, C - 4, fin_x);

		B++;

		if (B == B_MAX || mbgen + 1 == snpNum) { // process the block of SNPs using multi-threading
		
			vector <vector <double> > genoChrBlock (B, vector <double> (N)); //genoChrTrainingBlock.reserve(B_MAX);
			omp_set_num_threads(threads);
#pragma omp parallel for schedule(dynamic)
			for (int b = 0; b < B; b++) {

				int t = omp_get_thread_num();
				if (bufLens[b] > bufs[t].size()) bufs[t].resize(bufLens[b]);
				vector <double> genoChrVec;
				getSnpStatsBgen(&bufs[t][0], bufLens[b], &zBufs[b][0], zBufLens[b], lut,
					Nbgen, CVIndex, label, genoChrVec);
				genoChrBlock[b] = genoChrVec;
			}

			for (int i = 0; i < genoChrBlock.size(); i++) {
			
				genoChr.push_back(genoChrBlock[i]);
			}
			B = 0; // reset current block size
		}
	}
	free(allele0);
	free(allele1);
	return 0;
}

int STEPSEL::streamBgenSingle(vector <double> summInfoSingle, vector <int> CVIndex, string begin, string end, string label,
	vector <double> &genoVec)
{

	string bgen = begin + to_string((int)summInfoSingle[0]) + end;
	FILE* fin_x = fopen(bgen.c_str(), "rb");
	/********** READ HEADER **********/
	uint offset; fread(&offset, 4, 1, fin_x);
	uint L_H; fread(&L_H, 4, 1, fin_x);
	uint Mbgen; fread(&Mbgen, 4, 1, fin_x);
	uint Nbgen; fread(&Nbgen, 4, 1, fin_x);

	char magic[5]; fread(magic, 1, 4, fin_x); magic[4] = '\0';
	fseek(fin_x, L_H - 20, SEEK_CUR);
	uint flags; fread(&flags, 4, 1, fin_x);

	/********** READ SNP BLOCKS IN BATCHES **********/
	char snpID[65536], rsID[65536], chrStr[65536];
	char *allele1, *allele0;
	uint maxLA = 65536, maxLB = 65536;
	allele1 = (char *)malloc(maxLA + 1);
	allele0 = (char *)malloc(maxLB + 1);

	double lut[256];
	for (int i = 0; i <= 255; i++)
		lut[i] = i / 255.0;

	// during single-threaded reading of block, store SNP data for later multi-threaded processing
	//vector <uchar> bufs(1);
	vector <uchar> zBufs;
	uint zBufLens, bufLens;
	fseek(fin_x, offset + 4, SEEK_SET); //cout << "start: " << ftell(fin_x) << endl;
	fseek(fin_x, summInfoSingle[3], SEEK_CUR);//cout << "current: " << ftell(fin_x) << endl;
	ushort LS; fread(&LS, 2, 1, fin_x);
	fread(snpID, 1, LS, fin_x);
	ushort LR; fread(&LR, 2, 1, fin_x);
	fread(rsID, 1, LR, fin_x);
	ushort LC; fread(&LC, 2, 1, fin_x);
	fread(chrStr, 1, LC, fin_x); chrStr[LC] = '\0';
	uint physpos; fread(&physpos, 4, 1, fin_x); //cout << "physpos: " << physpos << endl;
	ushort K; fread(&K, 2, 1, fin_x); //cout << "K: " << K << endl;
	uint LA; fread(&LA, 4, 1, fin_x); //cout << "LA: " << LA << endl;
	fread(allele1, 1, LA, fin_x); //cout << "allele1: " << allele1 << endl;
	uint LB; fread(&LB, 4, 1, fin_x); //cout << "LB: " << LB << endl;
	fread(allele0, 1, LB, fin_x); //cout << "allele0: " << allele0 << endl;
	uint C; fread(&C, 4, 1, fin_x); //cout << "C: " << C << endl;
	if (C > zBufs.size()) zBufs.resize(C - 4);
	uint D; fread(&D, 4, 1, fin_x); //cout << "D: " << D << endl;
	zBufLens = C - 4; bufLens = D;
	fread(&zBufs[0], 1, C - 4, fin_x);
	// if (bufLens > bufs[0].size()) bufs[0].resize(bufLens);
	vector <uchar> bufs(bufLens);
	
	getSnpStatsBgen(&bufs[0], bufLens, &zBufs[0], zBufLens, lut, Nbgen, CVIndex, label, genoVec);
	
	free(allele0);
	free(allele1);
	return 0;
}

int STEPSEL::getGenotype1(vector<vector <string> > summData, float Len, float threshold1, float threshold2,
	int NTraining, int NTesting, vector<int> CVIndex, string begin, string end, int threads,
	string label, MatrixXd &genoTotMat, vector<vector <double> > &summInfoIdp)
{

	if (label != "Training" && label != "Testing") {

		cerr << "ERROR: The label is wrong!" << endl;
		exit(1);
	}
	int N;
	if (label == "Training") N = NTraining; 
	if (label == "Testing") N = NTesting; 

	//select snps in the threshold interval
	// vector<vector <double> > summInfo;
	vector <double> chrVec;
	vector<vector <double> > summInfo;
	for (int i = 0; i < summData.size(); i++) {

		if (atof(summData[i][2].c_str()) < threshold2 && atof(summData[i][2].c_str()) > threshold1) {

			vector <double> summFliter(6);
			summFliter[0] = atof(summData[i][0].substr(3).c_str());
			chrVec.push_back(summFliter[0]);
			summFliter[1] = atof(summData[i][1].c_str());
			summFliter[2] = atof(summData[i][2].c_str());
			summFliter[3] = atof(summData[i][3].c_str());
			summFliter[4] = atof(summData[i][4].c_str());
			summFliter[5] = atof(summData[i][5].c_str());
			summInfo.push_back(summFliter);
		}
	}
	cout << "The whole genome includes " << summInfo.size() << " snps < " << threshold2 << "." << endl;

	vector <double>::iterator chrVecIt;
	chrVecIt = unique(chrVec.begin(), chrVec.end());
	chrVec.resize(distance(chrVec.begin(), chrVecIt));

	int pos = 0;
	vector <vector <double> > genoTot;
	//chrVec.size()
	for (int i = 0; i < chrVec.size(); i++) {

		//get chromesome
		vector <vector <double> > summInfoChr, summInfoChrIdp;
		for (int j = pos; j < summInfo.size(); j++) {

			if (summInfo[j][0] == chrVec[i]) {

				summInfoChr.push_back(summInfo[j]);
			}
			else {

				pos = j+1;
				break;
			}
		}

		//get the information of independent snps
		snpSelChr(summInfoChr, Len, summInfoChrIdp);
		cout << "Chromosome " << chrVec[i] << " includes " << summInfoChr.size() << " snps < " << threshold2 << "." << endl;
		cout << "Chromosome " << chrVec[i] << " includes " << summInfoChrIdp.size() << " independent snps < " << threshold2 << "." << endl;

		//get the independent genotypes
		vector <vector <double> > genoChr;
		streamBgenChr(summInfoChrIdp, chrVec[i], N, CVIndex, begin, end, threads, label, genoChr);
		
		for (int i = 0; i < genoChr.size(); i++) {

			genoTot.push_back(genoChr[i]);
			summInfoIdp.push_back(summInfoChrIdp[i]);
		}
		cout << "Chromosome " << chrVec[i] << " is loaded!" << endl;
	}
	cout << "The whole genome includes " << genoTot.size() << " independent snps < " << threshold2 << "." << endl;
	//transformate vector to MatrixXd
	int varNum = genoTot.size();
	genoTotMat.resize(varNum, N);
	for (int i = 0; i < varNum; i++) {

		if (label == "Training") { genoTotMat.row(i) = standardizeVec(VectorXd::Map(&genoTot[i][0], N)); }
		if (label == "Testing") { genoTotMat.row(i) = standardizeVec(VectorXd::Map(&genoTot[i][0], N), summInfoIdp[i][4], summInfoIdp[i][5]); }
	}
	genoTotMat.transposeInPlace();

	return 0;
}

int STEPSEL::getGenotype2(vector<vector <string> > summData, float Len, float threshold1, float threshold2,
	int NTraining, int NTesting, vector<int> CVIndex, string begin, string end, int threads,
	string label, VectorXd resid, vector <double> &pTot, vector<vector <double> > &snpIdp)
{

	if (label != "Training" && label != "Testing") {

		cerr << "ERROR: The label is wrong!" << endl;
		exit(1);
	}
	int N;
	if (label == "Training") { N = NTraining; }
	if (label == "Testing") { N = NTesting; }

	//select snps in the threshold interval
	vector<vector <double> > summSel;
	vector <double> chrVec;
	for (int i = 0; i < summData.size(); i++) {

		if (atof(summData[i][2].c_str()) < threshold2 && atof(summData[i][2].c_str()) > threshold1) {

			vector <double> summFliter(6);
			summFliter[0] = atof(summData[i][0].substr(3).c_str());
			chrVec.push_back(summFliter[0]);
			summFliter[1] = atof(summData[i][1].c_str());
			summFliter[2] = atof(summData[i][2].c_str());
			summFliter[3] = atof(summData[i][3].c_str());
			summFliter[4] = atof(summData[i][3].c_str());
			summFliter[5] = atof(summData[i][3].c_str());
			summSel.push_back(summFliter);
		}
	}
	cout << "The whole genome includes " << summSel.size() << " snps < " << threshold2 << "." << endl;

	vector <double>::iterator chrVecIt;
	chrVecIt = unique(chrVec.begin(), chrVec.end());
	chrVec.resize(distance(chrVec.begin(), chrVecIt));

	int pos = 0;
	// chrVec.size()
	for (int i = 0; i < chrVec.size(); i++) {

		//get chromesome
		vector <vector <double> > snpChr;
		vector <vector <double> > snpChrIdp;
		for (int j = pos; j < summSel.size(); j++) {

			if (summSel[j][0] == chrVec[i]) {

				snpChr.push_back(summSel[j]);
			}
			else {

				pos = j + 1;
				break;
			}
		}

		//get the information of independent snps
		snpSelChr(snpChr, Len, snpChrIdp);
		cout << "Chromosome " << chrVec[i] << " includes " << snpChr.size() << " snps < " << threshold2 << "." << endl;
		cout << "Chromosome " << chrVec[i] << " includes " << snpChrIdp.size() << " independent snps < " << threshold2 << "." << endl;
		
		//get the independent genotypes
		vector <vector <double> > genoChr;
		streamBgenChr(snpChrIdp, chrVec[i], N, CVIndex, begin, end, threads, label, genoChr);
		cout << "Chromosome " << chrVec[i] << " is loaded!" << endl;
		MatrixXd genoChrMat(genoChr.size(), N);
		for (int i = 0; i < genoChr.size(); i++) {

			genoChrMat.row(i) = standardizeVec(VectorXd::Map(&genoChr[i][0], N));
		}
		genoChrMat.transposeInPlace();
		cout << "Transformation is ok!" << endl;
		vector <double> pChr;
		for (int i = 0; i < genoChrMat.cols(); i++) {

			double p = lmTest(genoChrMat.col(i), resid);
			pChr.push_back(p);
		}
		cout << "Single variable testing is ok!" << endl;
		for (int i = 0; i < pChr.size(); i++) {

			pTot.push_back(pChr[i]);
			snpIdp.push_back(snpChrIdp[i]);
		}
	}
	cout << "The whole genome includes " << pTot.size() << " independent snps < " << threshold2 << "." << endl;
	return 0;
}