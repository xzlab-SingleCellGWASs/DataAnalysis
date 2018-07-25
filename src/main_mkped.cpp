#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <sys/types.h>

#include "mkped.hpp"

using namespace std;

int main(int argc, char * argv[])
{
	MKPED cMkped;
	PARAM cPar;

	cMkped.Assign(argc, argv, cPar);
	cMkped.BatchRun(cPar);
	return EXIT_SUCCESS;
}