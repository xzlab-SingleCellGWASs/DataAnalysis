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


#ifndef __STEPPRED_H__
#define __STEPPRED_H__

#include <vector>
#include <string>
#include <iostream>

#include "Types.hpp"
#include "io.hpp"
#include "snpEff.hpp"
#include "mathFunc.hpp"

using namespace std;

class PARAM {

public:
	//parameters
	float prop;
	int seed;
	string phenoStr;
	string sqcNAStr;
	string summDataStr;
	double threshold1;
	double threshold2;
	double thresholdStep;
	float Len1;
	float Len2;
	int threads;	
	string bgenStart;
	string outpath;
	string snpfile;
	string phenofile;
};


class STEPRD {

public:
	//parameters
	string version;
	string date;
	string year;

	//constructor
	STEPRD(void);

	//functions
	void printHeader(void);
	void printHelp(void);
	void Assign(int argc, char ** argv, PARAM &cPar);
	void BatchRun(PARAM &cPar);
};

#endif