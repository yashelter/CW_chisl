#ifndef EULER_H
#define EULER_H

#include <vector>
#include <functional>
#include <utility>
#include <fstream>
#include "parameters.h"


class Euler
{
public:
	static std::vector<std::vector<double>> calculate_by_euler(const Parameters &parameters, double u0, double t0);
	static void print_results(std::ofstream& file, const std::vector<std::vector<double>>& data, const Parameters& parameters);
private:
	static std::vector<std::vector<double>> eulerCauchy(double u0, double t, double uf,
		std::function<double(double, double)> f,
		std::function<double(double)> fi);
};




#endif //EULER_H
