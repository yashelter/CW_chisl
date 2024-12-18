#include <euler.h>
#include <iostream>
#include "parameters.h"


int main()
{
	Globals globals;
	std::cout << globals.get_air_density(500) << std::endl;
	try
	{
		Parameters params = Parameters::load_from_file("/home/yashelter/Documents/CW_Numerical_Methods/src/resources/input.txt");

		double eps = 0.0001, h = 0.001;

		auto result = Euler::calculate_by_euler(params, params.u_0, params.t_0);
		std::ofstream outFile("result.txt");

		Euler::print_results(outFile, result, params);

	}
	catch (const std::exception& e)
	{
		std::cerr << e.what() << '\n';
		return 1;
	}


	return 0;
}
