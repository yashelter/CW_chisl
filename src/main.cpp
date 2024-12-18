#include <euler.h>
#include <iostream>
#include "parameters.h"


int main()
{
	try
	{
		Parameters params = Parameters::load_from_file("/home/yashelter/Documents/CW_chisl/src/resources/input.txt");
		auto result = Euler::calculate_by_euler(params, params.u_0, params.t_0);
		std::ofstream outFile("/home/yashelter/Documents/CW_chisl/src/result/data.txt");

		Euler::print_results(outFile, result, params);

	}
	catch (const std::exception& e)
	{
		std::cerr << e.what() << '\n';
		return 1;
	}


	return 0;
}
