#include <iostream>
#include <euler.h>

#include "parameters.h"

std::vector<std::vector<double>> Euler::eulerCauchy(double u0, double t0, double uf,
	std::function<double(double, double)> f,
	std::function<double(double)> fi)
{
	using std::abs;
	const double eps = 1e-5;
	double h = 1e-7;
	double uPrev = u0 + 2 * eps;
	double u = u0, t = t0, xi = 0; // Начальные значения
	// Вектор для хранения результатов (u, t, шаг h, расстояние xi, производная f(u, t))
	std::vector<std::vector<double>> result;

	result.push_back({u, t, 0, xi, fi(u)});
	int i = 0;

	while (std::abs(u - uPrev) > eps && abs(abs(uf) - abs(u)) > h) {
		i++;
		uPrev = u;
		double uPredict = u + h * f(u, t);
		double tPredict = t + h;

		long double val = (f(u, t) + f(uPredict, tPredict));
		// Итерационная обработка (коррекция)
		double uCorrect = u + (h / 2.0) * val;

		u = uCorrect;
		t = tPredict;

		xi = xi +  0.5 * h * (uPrev + u);
		result.push_back({u, t, h, xi, fi(u)});

		long double delta = std::abs((u - uPrev));
		if (delta > 1) {
			h *= 0.8;
		} else if (delta < 0.1) {
			h *= 1.2;
		}
	}

	return result;
}

std::vector<std::vector<double>> Euler::calculate_by_euler(const Parameters &parameters, double u0, double t0)
{
	std::function<double(double)> fi = [&parameters](double u)
	{
		return altha_u(parameters, u);
	};

    std::function<double(double, double)> f = [&parameters](double u, double t)
    {
        const double f_i = altha_u(parameters, u);
        return f_i * (parameters.u_f - u);
    };
    return eulerCauchy(u0, t0, parameters.u_f, f, fi);
}


void Euler::print_results(std::ofstream& file, const std::vector<std::vector<double>>& data, const Parameters& parameters)
{
	file << "# Things that we have\n";
	file << "# Final flow velocity: " << parameters.u_f << " (m/s) - The final velocity of the flow\n";
	file << "# Initial particle velocity: " << parameters.u_0 << " (m/s) - The initial velocity of the particle\n";
	file << "# Particle radius: " << parameters.r_is << " (m) - The radius of the particle\n";
	file << "# Nozzle radius: " << parameters.r_n << " (m) - The radius of the nozzle\n";
	file << "# Particle density: " << parameters.q << " (kg/m^3) - The density of the particle\n";
	file << "# Gas density: " << parameters.q_g << " (kg/m^3) - The density of the gas\n";
	file << "# Gas constant (R): " << parameters.R << " (J/(kg*K)) - The specific gas constant\n";
	file << "# M/M <- molar mass by molar mass: " << parameters.M_air << "of nothing\n";
	file << "# Molar mass of air (M): " << parameters.M << " (kg/mol) - The molar mass of air\n";
	file << "# Gas temperature (T): " << parameters.T << " (K) - The temperature of the gas\n";
	file << "# Reference temperature (T_0): " << parameters.T_0 << " (K) - The reference temperature for gas\n";
	file << "# Reference viscosity (eta_0): " << parameters.eta_0 << " (Pa*s) - The reference viscosity of the gas\n";
	file << "# Specific heat capacity of the gas (Cp): " << parameters.Cp << " (J/(kg*K)) - The specific heat capacity of the gas\n";
	file << "# Speed of sound (Sound speed): " << parameters.a << " (m/s) - The speed of sound in the gas\n";
	file << "# Dynamic (Viscosity): " << parameters.viscosity << " (Pa*s) - The dynamic viscosity of the gas\n";
	file << "\n\n\n";
	file << "# i  t  v  x  h  fi\n";
	//(u, t, шаг h, расстояние xi, производная f(u, t))

	long double distance = 0;

	for (int i = 0; i < data.size(); ++i)
	{
		distance += data[i][2];
		file << i <<  " " << data[i][1] << " " << data[i][0] << " " << data[i][3]  << " " << data[i][4] << "\n";
	}
	double max_t = data[data.size() - 1][1];

	file << "\n\n# Total time: " << max_t << " s";
	file << "\n# Total distance: " << distance  << " m";
	std::cout << "\nTotal time: " << max_t<< " s\n";

	file.close();
}
