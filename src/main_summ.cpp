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

#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <sys/types.h>

#include "summ.hpp"

using namespace std;

int main(int argc, char * argv[])
{
	SUMM cSumm;
	PARAM cPar;

	if (argc <= 1) {
		cSumm.printHeader();
		return EXIT_SUCCESS;
	}
	if (argc == 2 && argv[1][0] == '-' && argv[1][1] == 'h') {
		cSumm.printHelp();
		return EXIT_SUCCESS;
	}
	cSumm.Assign(argc, argv, cPar);
	cSumm.BatchRun(cPar);
	return EXIT_SUCCESS;
}