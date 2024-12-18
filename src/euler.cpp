#include <iostream>
#include <euler.h>

#include "parameters.h"

std::vector<std::vector<double>> Euler::eulerCauchy(double u0, double t0, std::function<double(double, double)> f) {
	const double eps = 1e-5;
	double h = 0.001;
	double uPrev = u0 + 2 * eps;
	double u = u0, t = t0, xi = 0; // Начальные значения
	// Вектор для хранения результатов (u, t, шаг h, расстояние xi, производная f(u, t))
	std::vector<std::vector<double>> result;
	f(u, t);
	result.push_back({u, t, 0, xi, f(u, t)});

	while (std::abs(u - uPrev) > eps) {
		uPrev = u;
		double uPredict = u + h * f(u, t);
		double tPredict = t + h;
		// Итерационная обработка (коррекция)
		double uCorrect = u + (h / 2.0) * (f(u, t) + f(uPredict, tPredict));

		u = uCorrect;
		t = tPredict;
		if (std::abs(f(u, t)) > 1) {
			h *= 0.8; // Уменьшаем шаг, если изменения слишком резкие
		} else if (std::abs(f(u, t)) < 0.1) {
			h *= 1.1; // Увеличиваем шаг, если изменения незначительны
		}
		xi += 0.5 * h * (uPrev + u);
		result.push_back({u, t, h, xi, f(u, t)});
	}

	return result;
}

std::vector<std::vector<double>> Euler::calculate_by_euler(const Parameters &parameters, double u0, double t0)
{
    std::function<double(double, double)> f_i = [&parameters](double u, double t)
    {
        const double f_i = altha_u(parameters, u);
        return f_i * (parameters.u_f - u);
    };
    return eulerCauchy(u0, t0, f_i);
}



void Euler::print_results(std::ofstream& file, const std::vector<std::vector<double>>& data, const Parameters& parameters)
{

	file << "# Constants and Parameters\n\n";
	file << "# Final flow velocity (u_f): " << parameters.u_f << " (m/s) - The final velocity of the flow\n";
	file << "# Initial particle velocity (u_0): " << parameters.u_0 << " (m/s) - The initial velocity of the particle\n";
	file << "# Particle radius (r_is): " << parameters.r_is << " (m) - The radius of the particle\n";
	file << "# Nozzle radius (r_n): " << parameters.r_n << " (m) - The radius of the nozzle\n";
	file << "# Particle density (q): " << parameters.q << " (kg/m^3) - The density of the particle\n";
	file << "# Gas density (q_g): " << parameters.q_g << " (kg/m^3) - The density of the gas\n";
	file << "# Gas constant (R): " << parameters.R << " (J/(kg*K)) - The specific gas constant\n";
	file << "# Gamma: " << parameters.gamma;
	file << "# Molar mass of air (M): " << parameters.M << " (kg/mol) - The molar mass of air\n";
	file << "# Gas temperature (T): " << parameters.T << " (K) - The temperature of the gas\n";
	file << "# Reference temperature (T_0): " << parameters.T_0 << " (K) - The reference temperature for gas\n";
	file << "# Reference viscosity (eta_0): " << parameters.eta_0 << " (Pa*s) - The reference viscosity of the gas\n";
	file << "# Specific heat capacity of the gas (Cp): " << parameters.Cp << " (J/(kg*K)) - The specific heat capacity of the gas\n";
	file << "# Speed of sound (Sound speed): " << parameters.a << " (m/s) - The speed of sound in the gas\n";
	file << "# Dynamic viscosity (Viscosity): " << parameters.viscosity << " (Pa*s) - The dynamic viscosity of the gas\n";
	file << "\n\n\n";


	file << "\n\n\n# i  t  v  x  h  fi\n";
	/*

	double max_t = data[data.size() - 1].first.first.first;
	double max_s = data[data.size() - 1].first.first.second;
	for (long i = 0; i < data.size(); ++i)
	{
		double fi =  data[i].second.second;
		double h =  data[i].second.first;
		double xi = data[i].first.first.second;
		double u = data[i].first.second.first;
		double t = data[i].first.second.second;
		file << i <<  " " << t << " " << u << " " << xi << " " << h << " " << fi <<"\n";
	}

	file << "\n\n# Total time: " << max_t << " s";
	file << "\n# Total distance: " << max_s  << " m";
	std::cout << "\nTotal time: " << max_t<< " s\n";
*/
	file.close();

}
