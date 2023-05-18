#include "ProductionPlan.h"

ProductionPlan::ProductionPlan(int periods, int blocks, int arcs, int units, int facilities)
{
	system("copy 05_debugger.txt predebugger.txt");
	remove("debugger.txt");

	debugFile.open("05_debugger.txt", ios::app);
	debugFile.setf(ios::fixed, ios::floatfield);
	debugFile.precision(2);

	nPeriods = periods;
	nBlocks = blocks;
	nArcs = arcs;
	uIndicator = units;
	nFacilities = facilities;

	metalPrice = refineCost = mineCost = discountRate = mineLowerCapacity = mineUpperCapacity = processLowerCapacity = processUpperCapacity = 0.0;

	mdouble.newMTX(bX, nBlocks);
	mdouble.newMTX(bY, nBlocks);
	mdouble.newMTX(bZ, nBlocks);
	mdouble.newMTX(bTons, nBlocks);
	mdouble.newMTX(bGrade, nBlocks);
	mdouble.newMTX(processCost, nFacilities); // Initialize processingCosts with the correct dimensions
	mdouble.newMTX(processRecovery, nBlocks, nFacilities); // Initialize processRecovery with the correct dimensions
	mdouble.newMTX(blockValue, nBlocks);

	mint.newMTX(blockID, nBlocks);

	mint.newMTX(fromBlock, nArcs);
	mint.newMTX(toBlock, nArcs);

}

ProductionPlan::~ProductionPlan(void)
{
	mdouble.delMTX(bX, nBlocks);
	mdouble.delMTX(bY, nBlocks);
	mdouble.delMTX(bZ, nBlocks);
	mdouble.delMTX(bTons, nBlocks);
	mdouble.delMTX(bGrade, nBlocks);
	mdouble.delMTX(blockValue, nBlocks);
	mdouble.delMTX(processCost, nFacilities); // Delete processingCosts with the correct dimensions
	mdouble.delMTX(processRecovery, nBlocks, nFacilities); // Delete processRecovery with the correct dimensions

	mint.delMTX(blockID, nBlocks);

	mint.delMTX(fromBlock, nArcs);
	mint.delMTX(toBlock, nArcs);
}

void ProductionPlan::ReadData()
{
	debugFile << "Reading XYZ and block grades ...."
		<< endl << endl << endl;

	std::cout << "Reading XYZ and block grades ....."
		<< endl << endl << endl;

	ifstream myInput("01_blockmodel.txt", ios::in);
	ifstream myPrecedence("02_precedence.txt", ios::in);

	if (myInput.is_open())
	{
		for (int i = 0; i < nBlocks; i++)
		{
			myInput >> bX[i] >> bY[i] >> bZ[i] >> bGrade[i] >> bTons[i];
			for (int p = 0; p < nFacilities; p++) // Read the processing recovery for each block and each processing facility
            {
                myInput >> processRecovery[i][p]; // where i is the block and p is the corresponding Processing Facility
            }
			//std::cout << (i + 1) << endl;
		}
	}

	myInput.close();

	if (myPrecedence.is_open())
	{
		for (int i = 0; i < nArcs; i++)
		{
			myPrecedence >> fromBlock[i] >> toBlock[i];
			//std::cout << (i + 1) << endl;
		}
	}

	myPrecedence.close();

	debugFile << "Reading XYZ and block grades ...."
		<< endl << endl << endl;

	std::cout << "Reading XYZ and block grades ....."
		<< endl << endl << endl;
}

void ProductionPlan::AssignBlockIDS()
{
	debugFile << "Assign IDs ..."
		<< endl << endl << endl;

	std::cout << "Assign IDs ....."
		<< endl << endl << endl;

	int k = 0;

	for (int j = 0; j < nBlocks; j++)
	{
		blockID[k] = k + 1;
		k++;
	}

	ofstream myOutput("03_InitialData.txt", ios::out);
	myOutput.setf(ios::fixed, ios::floatfield);

	    for (int i = 0; i < nBlocks; i++)
    {
        myOutput << setiosflags(ios::left) << setprecision(2)
            << setw(10) << blockID[i]
            << setw(10) << bX[i]
            << setw(10) << bY[i]
            << setw(10) << bZ[i]
            << setw(10) << bTons[i]
            << setw(10) << bGrade[i];
        
        for (int p = 0; p < nFacilities; p++)
        {
            myOutput << setw(10) << processRecovery[i][p];
        }

        myOutput << endl;
    }

    for (int i = 0; i < nArcs; i++)
    {
        myOutput << setiosflags(ios::left)
            << setw(10) << fromBlock[i]
            << setw(10) << toBlock[i] << endl;
    }

	debugFile << "Assign IDs ...."
		<< endl << endl << endl;

	std::cout << "Assign IDs ....."
		<< endl << endl << endl;
}

void ProductionPlan::AssignData(double mP, double rC, double mC, double dR, double mLC, double mUC, double pLC, double pUC, double* pC)
{
	metalPrice = mP;
	refineCost = rC;
	mineCost = mC;
	discountRate = dR;
	mineLowerCapacity = mLC;
	mineUpperCapacity = mUC;
	processLowerCapacity = pLC;
	processUpperCapacity = pUC;
	//assign processing costs
	for (int p = 0; p < nFacilities; p++)
		{
		processCost[p] = pC[p];
		}
}

void ProductionPlan::CalculateBlockValues() 
{
	for (int i = 0; i < nBlocks; i++)
	{
		double maxValue = -std::numeric_limits<double>::max();

		for (int p = 0; p < nFacilities; p++)
		{
			double value = 0.0;
			if (uIndicator == 1)
			{
				value = ((metalPrice - refineCost) * (bGrade[i] / 100) * processRecovery[i][p] * bTons[i]) - (mineCost * bTons[i]) - (processCost[p] * bTons[i]);
			}
		
			else if (uIndicator == 2)
			{
				value = ((metalPrice - refineCost) * bGrade[i] * processRecovery[i][p] * bTons[i]) - (mineCost * bTons[i]) - (processCost[p] * bTons[i]);
			}
			if (value > maxValue)
			{
				maxValue = value;
			}
		}
			if (maxValue <= 0.0)
			{
				maxValue = (-1) * mineCost * bTons[i];
			}

			blockValue[i] = maxValue;
		
		debugFile << setiosflags(ios::left) << setw(10) << (i + 1) << blockValue[i] << endl;
	}

	debugFile << endl;
}

void ProductionPlan::CreateModel()
{
	IloEnv env;
	IloModel mod(env);
	IloTimer solutiontime(env);

	typedef IloArray<IloArray<IloNumVarArray> > Array3D;

	Array3D x(env, nBlocks); // Decision variable for the production of each block at each facility in each period

	for(int i = 0; i < nBlocks; i++) 
	{
		x[i] = IloArray<IloNumVarArray>(env, nFacilities); 

		for (int p = 0; p < nFacilities; p++)
		{

		x[i][p] = IloNumVarArray(env, nPeriods, 0, 1, ILOINT); 

		}
	}
	
	for (int i = 0; i < nBlocks; i++)
	{
		for (int p = 0; p < nFacilities; p++)
		{
			for (int t = 0; t < nPeriods; t++) 

			debugFile << i << "  " << p << " " << t << "  " << x[i][p][t] << endl; 
		}
	}

	debugFile << endl << endl;

	IloExpr objExpression(env); //Objective Function Expression

	int n = 0; //Loop integer, Loops until the objective completed.

	for (int i = 0; i < nBlocks; i++)
	{
		for (int p = 0; p < nFacilities; p++)
		{
			for (int t = 0; t < nPeriods; t++)
			{

			n++;

			objExpression += (blockValue[i]/(pow((1+(discountRate/100)),(t+1)))) * x[i][p][t]; //Objective Function Maximising the NPV        

			debugFile << n << "  " << i << "  " << p << "  " << t << "  " << x[i][p][t] << endl;
			}
		}		
	}

	debugFile << endl << endl;

	IloObjective objFunction(env, objExpression, IloObjective::Maximize);

	mod.add(objFunction);
	
	for (int i = 0; i < nBlocks; i++) 
	{
		IloExpr reserveExpression(env);

		for (int t = 0; t < nPeriods; t++) 
		{
			for (int p = 0; p < nFacilities; p++) /*Should I loop around facilities? Would it be a problem?*/
			{
				reserveExpression += x[i][p][t];
			}
		}

		mod.add(reserveExpression <= 1); // Reserve Constraint - Each Block can be mined at once!
	}
	 /*//Processing Facility Selection Constraint - Each block can assign only one processing facility
	for (int i = 0; i < nBlocks; i++)
	{
		for (int t = 0; t < nPeriods; t++)
		{
			IloExpr facilitySelection(env);

			for (int p = 0; p < nFacilities; p++)
			{
				facilitySelection += x[i][p][t];
			}

			mod.add(facilitySelection == 1);
		}
	}*/
	

	for (int t = 0; t < nPeriods; t++)
{
    IloExpr miningExpression(env);

    for (int i = 0; i < nBlocks; i++)
    {
        for (int p = 0; p < nFacilities; p++) 
        {
            miningExpression += (bTons[i] * x[i][p][t]); 
        }
    }

    mod.add(miningExpression >= mineLowerCapacity); // Mining Capacity Constraint Lower Limit
    mod.add(miningExpression <= mineUpperCapacity); 
}
	
	for (int t = 0; t < nPeriods; t++)
	{
		for (int p = 0; p < nFacilities; p++)
		{
			IloExpr processExpression(env);

        for (int i = 0; i < nBlocks; i++)
        {
            if (blockValue[i] > 0)
            {
                processExpression += (bTons[i] * x[i][p][t]); // Processing Capacities should differ for each facility!!!!!
            }
        }

        mod.add(processExpression >= processLowerCapacity); // Lower Processing Capacity Constraint
		mod.add(processExpression <= processUpperCapacity); // Upper Processing Capacity Constraint
		}
	}
	   
	int fBlock = 0, tBlock = 0;

	for (int j = 0; j < nArcs; j++)
	{
		fBlock = fromBlock[j];
		tBlock = toBlock[j];

		debugFile << "fBlock = " << fBlock << "  " << fBlock - 1 << "   tBlock = " << tBlock << "  " << tBlock - 1 << endl;
 	 
		IloExpr precedenceTo(env);

		for (int t = 0; t < nPeriods; t++)
		{
			debugFile << x[fBlock - 1][t] << "  " << x[tBlock - 1][t] << endl;

			precedenceTo += x[tBlock - 1][t]; //Precedence Constraint 
			mod.add(x[fBlock - 1][t] - precedenceTo <= 0);

		debugFile << endl;
		}

	}

	IloCplex cplex(env);
	cplex.extract(mod);
	cplex.exportModel("04_Formulation.lp");
	
	cplex.setParam(IloCplex::ClockType, 1);
	cplex.setParam(IloCplex::EpGap, 0.01);
	cplex.setParam(IloCplex::TiLim, 2592000);

	auto start = std::chrono::system_clock::now();
	solutiontime.start();
	cplex.solve();
	solutiontime.stop();
	auto end = std::chrono::system_clock::now();

	std::chrono::duration<double> elapsed_seconds = end - start;

	ofstream myOutput("05_OptimalPP.txt", ios::out);
	myOutput.setf(ios::fixed, ios::floatfield);

	if (myOutput.is_open())
	{
		IloNum objValue = cplex.getObjValue();

		try {
			objValue = cplex.getObjValue();
		}
		catch (IloCplex::Exception& e) {
			std::cerr << "CPLEX Exception: " << e.getMessage() << std::endl;
			// handle the exception or exit the program
		}

		myOutput << "Discounted value of the production plan = " << setprecision(0) << objValue << endl;
		myOutput << "Solution time = " << setprecision(5) << solutiontime.getTime() << " seconds" << endl << endl;
		myOutput << "Solution time = " << setprecision(5) << elapsed_seconds.count() << " seconds" << endl << endl;

		for (int t = 0; t < nPeriods; t++)
		{
			for (int i = 0; i < nBlocks; i++)
			{
				for (int p = 0; p < nFacilities; p++) // Added this loop for processing facilities

					if (cplex.getValue(x[i][p][t]) > 0)
					{
						myOutput << setiosflags(ios::left) << setprecision(0)
							<< setw(10) << bX[i]
							<< setw(10) << bY[i]
							<< setw(10) << bZ[i]
							<< setw(10) << bTons[i]
							<< setprecision(5)
							<< setw(15) << bGrade[i]
							<< setw(10) << setprecision(0) << (t+1)
							<< setw(10) << p + 1 // Added facility number (p + 1) to the output
							<< setw(10) << cplex.getValue(x[i][p][t]) << endl;
				}
			}
		}

		ofstream myOutput1("06_PitProcess.txt", ios::out);
		myOutput1.setf(ios::fixed, ios::floatfield);

		for (int i = 0; i < nBlocks; i++)
		{
			for (int p = 0; p < nFacilities; p++) // Added this loop for processing facilities
			{
				for (int t = 0; t < nPeriods; t++)
				{
					if (cplex.getValue(x[i][p][t]) > 0)
					{
						myOutput1 << setiosflags(ios::left) << setprecision(0)
							<< setw(10) << blockID[i]
							<< setw(10) << (t) << endl;
					}
				}
			}			
		}
	}
}

void ProductionPlan::AnalysisOfResults()
{
	ofstream resultsFile("07_results.txt", ios::out);
	resultsFile.setf(ios::fixed, ios::floatfield);

	resultsFile << setiosflags(ios::left) << setprecision(0);

	ofstream modelFile("08_blockmodelSchedule.txt", ios::out);
	modelFile.setf(ios::fixed, ios::floatfield);

	ofstream blockToFacilityFile("09_blockToFacility.txt", ios::out);
	blockToFacilityFile.setf(ios::fixed, ios::floatfield);

	ofstream processedMaterialFile("10_processedMaterial.txt", ios::out);
	processedMaterialFile.setf(ios::fixed, ios::floatfield);

	//Pit to Process
	ifstream pToP("06_PitProcess.txt", ios::in);

	int *mbID, *mTime, *mFacility;

	mint.newMTX(mbID, (nBlocks * nPeriods));
	mint.newMTX(mTime, (nBlocks * nPeriods));
	mint.newMTX(mFacility, (nBlocks * nFacilities)); 
	int i = 0, nPitProcess = 0;
	while (pToP >> mbID[i] >> mTime[i]>> mFacility[i])
	{
		i++;
	}

	nPitProcess = i;

	for (int b = 0; b < nBlocks; b++)
	{
		for (int i = 0; i < nPitProcess; i++)
		{
			if (blockID[b] == mbID[i])
			{
				modelFile << setiosflags(ios::left) << setprecision(0)
					<< setw(10) << bX[b]
					<< setw(10) << bY[b]
					<< setw(10) << bZ[b]
					<< setw(10) << bTons[b]
					<< setprecision(5)
					<< setw(15) << bGrade[b]
					<< setw(15) << (mTime[i] + 1) << endl
					<< setw(15) << (mFacility[i] + 1) << endl;

				goto out_block;
			}
		}
	
		modelFile	<< setiosflags(ios::left) << setprecision(0)
					<< setw(10) << bX[b]
					<< setw(10) << bY[b]
					<< setw(10) << bZ[b]
					<< setw(10) << bTons[b]
					<< setprecision(5)
					<< setw(15) << bGrade[b]
					<< setw(15) << (nPeriods + 100) << endl;

		out_block: debugFile << endl;
	}

	double *qtWaste, *qtOre, *qtMetal, *cutG, *avgG, *cashF, *netPV, *totalProcessed;

	mdouble.newMTX(qtWaste, nPeriods);
	mdouble.newMTX(qtOre, nPeriods);
	mdouble.newMTX(qtMetal, nPeriods);
	mdouble.newMTX(cutG, nPeriods);
	mdouble.newMTX(avgG, nPeriods);
	mdouble.newMTX(cashF, nPeriods);
	mdouble.newMTX(netPV, nPeriods);
	mdouble.newMTX(totalProcessed, nFacilities);

	for (int p = 0; p < nFacilities; p++)
	{
		totalProcessed[p] = 0.0;
	}

	for (int i = 0; i < nPitProcess; i++) // loop over pit process entries
	{
		int b = mbID[i];
		int p = mFacility[i];
		int t = mTime[i];

		totalProcessed[p] += bTons[b];

		blockToFacilityFile << setiosflags(ios::left) << setprecision(0)
			<< setw(10) << "Block " << blockID[b] << "is sent to Processing Facility "
			<< setw(10) << (p + 1)<< "in period "
			<< setw(10) << (t + 1) << "." << endl;
	}

	for (int p = 0; p < nFacilities; p++)
	{
		processedMaterialFile << setiosflags(ios::left) << setprecision(0)
			<< setw(10) << "Processing Facility " << (p + 1) << " processes "
			<< setw(10) << totalProcessed[p] << " tons of material" << endl;
	}

	for (int t = 0; t < nPeriods; t++)
	{
			qtWaste[t] = 0.0;
			qtOre[t] = 0.0;
			qtMetal[t] = 0.0;
			cutG[t] = 1000000000.0;
			avgG[t] = 0.0;	
			netPV[t] = 0.0;
	}

	for (int i = 0; i < nPitProcess; i++) // Instead of the below code looping over all blocks,periods,facilities; loop over pit process entries
	{
		double sumG = 0.0; 
		int b = mbID[i];
		int p = mFacility[i];
		int t = mTime[i];

		if (blockValue[b] <= 0)
		{
			qtWaste[t] += bTons[b];
		}
		else
		{
			qtOre[t] += bTons[b];

			for (int p = 0; p < nFacilities; p++)
			{
				qtMetal[t] += (processRecovery[b][p] * bTons[b]);
			}
			sumG += (bGrade[b] * bTons[b]);

			if (bGrade[b] < cutG[t])
			{
				cutG[t] = bGrade[b];
			}
		}
	}
	/*for (int t = 0; t < nPeriods; t++)
	{
		double sumG = 0.0;

		for (int b = 0; b < nBlocks; b++)
		{
			for (int i = 0; i < nPitProcess; i++)
			{
				if (blockID[b] == mbID[i] && t == mTime[i])
				{
					if (blockValue[b] <= 0)
					{
						qtWaste[t] += bTons[b];
					}
					else
					{	
						qtOre[t] += bTons[b];

						for (int p = 0; p < nFacilities; p++)
						{
							qtMetal[t] += (processRecovery[b][p] * bTons[b]);
						}
						sumG += (bGrade[b] * bTons[b]);

						if (bGrade[b] < cutG[t])
						{
							cutG[t] = bGrade[b];
						}
					}
				}
			}
		}*/
	for (int t = 0; t < nPeriods; t++)
	{
		double sumG = 0.0;
		avgG[t] = (sumG / qtOre[t]);

		if (uIndicator == 1)
		{
			qtMetal[t] = (qtMetal[t] * (avgG[t] / 100));
		}

		if (uIndicator == 2)
		{
			qtMetal[t] = (qtMetal[t] * avgG[t]);
		}

		double totalProcessCost = 0.0;
		for (int i = 0; i < nPitProcess; i++) // loop over pit process entries
		{
			int b = mbID[i];
			int p = mFacility[i];

			if(t == mTime[i])
			{
				totalProcessCost += (processCost[p] * bTons[b]);
			}
		}

		cashF[t] = ((metalPrice - refineCost) * qtMetal[t]) - mineCost * (qtWaste[t] + qtOre[t]) - totalProcessCost;
	}

	for (int t = 0; t < nPeriods; t++)
	{
		for (int t1 = t; t1 < nPeriods; t1++)
		{
			if (t == 0)
			{
				netPV[t] += (cashF[t1] / (pow((1 + (discountRate / 100)), (t1 + 1))));
			}
			else
			{
				netPV[t] += (cashF[t1] / (pow((1 + (discountRate / 100)), (t1 - t + 1))));
			}			
		}
	}

	resultsFile << endl;

	resultsFile << setw(10) << "Year" 
			    << setw(20) << "Quantity of Waste" 
				<< setw(20) << "Cut-off Grade"
				<< setw(20) << "Average Grade"
				<< setw(20) << "Quantity of Ore" 
				<< setw(20) << "Quantity of Metal" 
				<< setw(20) << "Cash Flow" 
				<< setw(20) << "NPV" << endl;

	resultsFile << setw(10) << ""
				<< setw(20) << "(tonnes)"
				<< setw(20) << "(%)"
				<< setw(20) << "(%)"
				<< setw(20) << "(tonnes)" 
				<< setw(20) << "(tonnes)" 
				<< setw(20) << "($)" 
				<< setw(20) << "($)" << endl << endl;

	for (int t = 0; t < nPeriods; t++)
	{
		resultsFile << setw(10) << (t+1)
					<< setw(20) << qtWaste[t]
					<< setw(20) << setprecision(5) << cutG[t]
					<< setw(20) << setprecision(5) << avgG[t]
					<< setw(20) << setprecision(0) << qtOre[t] 
					<< setw(20) << setprecision(0) << qtMetal[t] 
					<< setw(20) << setprecision(0) << cashF[t] 
					<< setw(20) << setprecision(0) << netPV[t] << endl;
	}

	mint.delMTX(mbID, (nBlocks * nPeriods));
	mint.delMTX(mTime, (nBlocks * nPeriods));
	mint.delMTX(mFacility, (nBlocks * nPeriods));

	mdouble.delMTX(qtWaste, nPeriods);
	mdouble.delMTX(qtOre, nPeriods);
	mdouble.delMTX(qtMetal, nPeriods);
	mdouble.delMTX(cutG, nPeriods);
	mdouble.delMTX(avgG, nPeriods);
	mdouble.delMTX(cashF, nPeriods);
	mdouble.delMTX(totalProcessed, nFacilities);
	
}
