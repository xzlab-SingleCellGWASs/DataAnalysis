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


#ifndef __IO_H__
#define __IO_H__

#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <boost/random.hpp>

#include "Types.hpp"

using namespace std;

int readTable(string infile, int skip, char *separator, vector< vector<string> > &data);
int baseToCode(string base);
string codeToBase(int code);
string standardOutputSummaryData(vector <double> snpInfo, string snpName);
string standardOutputInfo(vector <double> snpInfo);
int getCVIndex(float prop, int seed, vector<vector <string> > sqcNAIndex, vector<vector <string> > pheno, vector <int> &CVIndex);

#endif