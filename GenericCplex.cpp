#include "ProductionPlan.h"

//remove "#include <ilcplex/ilocplex.h>" from the fist line and add to new class.

using namespace std; 

int main(){ 
 
	std::cout << "In main" << std::endl;

	ifstream file("00_baseParameters.txt", ios::in);
	ofstream outfile("progTime.txt", ios::out);

	int periods = 0; // Number of periods
	int blocks = 0; // Number of blocks
	int arcs = 0; // Number of precedence arcs
	int units = 0; // Number of mining units
	int facilities = 0; // Number of processing facilities
	double metal_price = 0.0, refining_cost = 0.0, mining_cost = 0.0, discount_rate = 0.0; 
	double mining_cap_low = 0.0, mining_cap_up = 0.0, processing_cap_low = 0.0, processing_cap_up = 0.0;

	file >> periods; // Read the number of periods from the input file
	file >> blocks; // Read the number of blocks from the input file
	file >> arcs; // Read the number of precedence arcs from the input file
	file >> units; // Read the number of mining units from the input file
	file >> facilities; // Read the number of processing facilities from the input file
		double* processCost = new double[facilities]; // Create an array for processing costs
	file >> metal_price >> refining_cost >> mining_cost >> discount_rate; // Read the metal price, refining cost, mining cost and discount rate from the input file
	file >> mining_cap_low >> mining_cap_up; // Read the mining capacity lower and upper bounds from the input file
	file >> processing_cap_low >> processing_cap_up; // Read the processing capacity lower and upper bounds from the input file
	
    for (int p = 0; p < facilities; p++) 
	{ 
        file >> processCost[p]; // Read processing costs for each facility from the input file
    }

	outfile << "Number of Blocks: " << blocks << endl; 
	outfile << "Number of precedence arcs: " << arcs << endl;

	ProductionPlan pp(periods, blocks, arcs, units, facilities); // Update the constructor with the number of processing facilities

	pp.ReadData();
	pp.AssignBlockIDS();

	pp.AssignData(metal_price, refining_cost, mining_cost, discount_rate, mining_cap_low, mining_cap_up, processing_cap_low, processing_cap_up, processCost);

	pp.CalculateBlockValues();

	pp.CreateModel();
	pp.AnalysisOfResults();

	delete[] processCost; // Delete the processing cost array

return 0;
}
