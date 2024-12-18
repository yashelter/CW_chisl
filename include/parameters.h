#ifndef CW_NM_PARAMETERS_H
#define CW_NM_PARAMETERS_H

#include <string>
#include <fstream>
#include <cmath>

struct Parameters
{

	const double P = 101325;
	const double R =  8.31446261815324;

	const double u_f;					 // Скорость потока
	const double u_0;					 // Начальная скорость частицы
	const double t_0 = 0.0;				 // Начальное время
	const double r_is;					 // Радиус частицы
	const double r_n;					 // Радиус сопла
	const double q = 3990;				 // Плотность частицы
	double q_g;							 // Плотность газа
	const double M_air = 1.4;		     // Отношение молярных теплоемкостей
	const double M = 0.029;		  	     // Молярная масса воздуха
	const double T;						 // Температура газа
	const double T_0 = 291.15;			 // Контрольная температура газа
	const double eta_0 = 1.827e-5;		 // Контрольная вязкость газа
	const double Cp = 920;				 // Удельная теплоемкость вещества

	double a; // Скорость звука
	double viscosity;    // Вязкость

	static Parameters load_from_file(const std::string& filename);

	Parameters(double u_f, double u_0, double r_is, double r_n, double T)
			: u_f(u_f), u_0(u_0), r_is(r_is), r_n(r_n), T(T)
	{
		q_g = calculate_gas_desity();
		a = (calculate_sound_speed(R, M_air, M, T));
		viscosity = (calculate_viscosity(eta_0, T_0, T));
	}

private:
	double calculate_gas_desity()
	{
		return P * M / (R * T);
	}

	static double calculate_sound_speed(double R, double gamma, double M, double T)
	{
		return std::sqrt((gamma * R * T) / M);
	}

	static double calculate_viscosity(double eta_0, double T_0, double T)
	{
		return eta_0 * std::pow(T / T_0, 1.5);
	}
};

double altha_u(const Parameters& params, const double& u);

double calculate_stocks (double q, double r_is, double a_s, double eta, double r_n);

double calculate_C_i (double M, double Re, double S, double u, double a, const Parameters& p);

double calculate_Mach_number (double u_f, double u, double a_s);

double calculate_Reynolds_number (double q, double r_is, double u, double u_f, double eta);

double calculate_S (double M, double gamma);


#endif //CW_NM_PARAMETERS_H
