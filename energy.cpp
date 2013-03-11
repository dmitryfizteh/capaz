#include "defines.h"

// !!! Нужно потом будет вынести в структуру констант
double T_0 = 273.; // К
double ro_r = 2000.; // кг/м^3
double lambda0_w = 0.553; // Вт/(м*К)
double lambda0_n = 0.14;
double lambda0_g = 0.0237;
double lambda0_r = 1.;
double lambdaA_w = 3E-3; // 1/K
double lambdaA_n = 1E-3;
double lambdaA_g = 0.82;
double c0_w = 4.194E3; // Дж/(кг*К)
double c0_n = 1.7E3;
double c0_g = 1E3;
double c0_r = 0.8E3;
double C_w = 1.15;
double C_w2 = 0.015;
double C_n = 3.4;
double C_g = 0.119;
double C_r = 0.75;

// Коэффициенты удельных теплоемкостей при постоянном давлении  для water, napl, gas and rock в Вт/(м*К)
double с_w (double T, consts def)
{
	return c0_w - C_w * (T - T_0) + C_w2 * (T - T_0) * (T - T_0);
}

double с_n (double T, consts def)
{
	return c0_n + C_n * (T - T_0);
}

double с_g (double T, consts def)
{
	return c0_g + C_g * (T - T_0);
}

double с_r (double T, consts def)
{
	return c0_r + C_r * (T - T_0);
}

// Коэффициенты теплопроводности для water, napl, gas and rock
double lambda_w (double T, consts def)
{
	return lambda0_w * (1 - lambdaA_w * (T - T_0));
}

double lambda_n (double T, consts def)
{
	return lambda0_n * (1 - lambdaA_n * (T - T_0));
}

double lambda_g (double T, consts def)
{
	return lambda0_g * pow((T / T_0), lambdaA_g);
}

double lambda_r (double T, consts def)
{
	return lambda0_r;
}

// Эффективный коэффициент теплопроводности в точке (будет использоваться при расчете теплового потока)
double assign_lambda_eff (ptr_Arrays HostArraysPtr, int local, consts def)
{
	return HostArraysPtr.m[local] * (HostArraysPtr.S_w[local] * lambda_w (HostArraysPtr.T[local], def)
		+ HostArraysPtr.S_n[local] * lambda_n (HostArraysPtr.T[local], def)
		+ HostArraysPtr.S_g[local] * lambda_g (HostArraysPtr.T[local], def)) 
		+ (1. - HostArraysPtr.m[local]) * lambda_r (HostArraysPtr.T[local], def);
}

// Расчет энтальпии по температуре и теплоемкости
// !!! Переписать, задав точность для метода Симпсона и передавая указатель на функцию, чтобы не копировать одно и то же
double assign_H_w (double T, consts def)
{
	/* Возможно, что Симпсон понадобиться позже, а пока все равно нужно знать явную зависимость энергии от температуры
	double integral = 0, sum = 0, h_temp;
	int N_temp = 50;

	h_temp = (T - T_0) / N_temp;
	

	integral += с_w(T_0, def);
	integral += с_w(T, def);

	for(int i = 2; i < N_temp; i+=2)
		sum += с_w(T_0 + i * h_temp, def);

	sum *= 2;
	integral += sum;
	sum = 0;

	for(int i = 1; i < N_temp; i+=2)
		sum += с_w(T_0 + i * h_temp, def);

	sum *= 4;
	integral += sum;

	h_temp /= 3;
	integral *= h_temp;

	return integral;
	*/
	return (T - T_0) * (c0_w - C_w * (T - T_0) / 2 + C_w2 * (T - T_0) * (T - T_0) / 3);
}

double assign_H_n (double T, consts def)
{
	return (T - T_0) * (c0_n + C_n * (T - T_0) / 2);
}

double assign_H_g (double T, consts def)
{
	return (T - T_0) * (c0_g + C_g * (T - T_0) / 2);
}

double assign_H_r (double T, consts def)
{
	return (T - T_0) * (c0_r + C_r * (T - T_0) / 2);
}

void assign_H (ptr_Arrays HostArraysPtr, int local, consts def)
{
	HostArraysPtr.H_w[local] = assign_H_w (HostArraysPtr.T[local], def);
	HostArraysPtr.H_n[local] = assign_H_n (HostArraysPtr.T[local], def);
	HostArraysPtr.H_g[local] = assign_H_g (HostArraysPtr.T[local], def);
	HostArraysPtr.H_r[local] = assign_H_r (HostArraysPtr.T[local], def);

	test_positive(HostArraysPtr.H_w[local], __FILE__, __LINE__);
	test_positive(HostArraysPtr.H_n[local], __FILE__, __LINE__);
	test_positive(HostArraysPtr.H_g[local], __FILE__, __LINE__);
	test_positive(HostArraysPtr.H_r[local], __FILE__, __LINE__);
}

// Коэффициенты вязкости для water, napl, gas and rock

// Расчет потока энергии в точке

// Расчет внутренней энергии всей системы в точке
void assign_E_current (ptr_Arrays HostArraysPtr, int local, consts def)
{
	HostArraysPtr.E[local] = (HostArraysPtr.m[local] * (HostArraysPtr.S_w[local] * (HostArraysPtr.ro_w[local] * HostArraysPtr.H_w[local] - HostArraysPtr.P_w[local])
		+ HostArraysPtr.S_n[local] * (HostArraysPtr.ro_n[local] * HostArraysPtr.H_n[local] - HostArraysPtr.P_n[local])
		+ HostArraysPtr.S_g[local] * (HostArraysPtr.ro_g[local] * HostArraysPtr.H_g[local] - HostArraysPtr.P_g[local])) 
		+ (1. - HostArraysPtr.m[local]) * (ro_r * HostArraysPtr.H_r[local] - HostArraysPtr.P_w[local]));

	test_positive(HostArraysPtr.E[local], __FILE__, __LINE__);
}

// Расчет теплового потока в точке
