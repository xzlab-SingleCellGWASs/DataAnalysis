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


#ifndef __SUMM_H__
#define __SUMM_H__

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
	float bgenMinMAF;
	float bgenMinINFO;
	float bgenMinHWE;
	float bgenMinCallingRate;
	int threads;
	float prop;
	int seed;
	string phenoStr;
	string sqcNAStr;
	string bgenStart;
	string outpath;
	string outfile;
	string chr;
	string cv;
};


class SUMM {

public:
	//parameters
	string version;
	string date;
	string year;

	//constructor
	SUMM(void);
	
	//functions
	void printHeader(void);
	void printHelp(void);
	void Assign(int argc, char ** argv, PARAM &cPar);
	void BatchRun(PARAM &cPar);
};

#endif
