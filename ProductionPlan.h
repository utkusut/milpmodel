#pragma once

#include "matrix_t.h"

#include <ilcplex/ilocplex.h>
#include <string> 
#include <iostream>
#include <ios>
#include <iomanip>
#include <fstream> 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sstream>
#include <vector>
#include <algorithm>
#include <time.h>
#include <ctime>
#include <chrono>
#include <ctype.h>
#include <assert.h>
#include <limits>

class ProductionPlan
{
public:
	ProductionPlan(int, int, int, int, int); // Added integer for the Number of Facilities
	~ProductionPlan();

	void ReadData();
	void AssignBlockIDS();
	void AssignData(double, double, double, double, double, double, double, double, double*); //Updated according to accept the vector format of Processing Costs
	void CalculateBlockValues();
	void CreateModel();
	void AnalysisOfResults();

protected:

	MATRIX<int> mint;
	MATRIX<double> mdouble;
	MATRIX<char> mchar;

	double *bX, *bY, *bZ, *bTons, *bGrade, *blockValue, *processCost;
	int *blockID, *fromBlock, *toBlock;
	double **processRecovery; // Change processRecovery to a two-dimensional vector- pointer to pointer

	char *debugfilename;
	ofstream debugFile;

private:
	int nPeriods, nBlocks, nArcs, uIndicator, nFacilities;

	double metalPrice, refineCost, mineCost, discountRate, mineLowerCapacity, mineUpperCapacity, processLowerCapacity, processUpperCapacity;
};

