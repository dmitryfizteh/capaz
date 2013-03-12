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

// Расчет теплового потока в точке
double assign_T_flow (ptr_Arrays HostArraysPtr, int i, int j, int k, consts def)
{	
	if ((i != 0) && (i != (def.locNx) - 1) && (j != 0) && (j != (def.locNy) - 1) && (((k != 0) && (k != (def.locNz) - 1)) || ((def.locNz) < 2)))
	{
		double T_flow = 0;
		int local=i + j * (def.locNx) + k * (def.locNx) * (def.locNy);

		if ((def.locNz) > 2)
		{
			T_flow += (assign_lambda_eff(HostArraysPtr, local + 1, def) * HostArraysPtr.T[local + 1]
			- 2 * assign_lambda_eff(HostArraysPtr, local, def) * HostArraysPtr.T[local]
			+ assign_lambda_eff(HostArraysPtr, local - 1, def) * HostArraysPtr.T[local - 1]) / ((def.hx) * (def.hx));
		}
		if ((def.locNz) > 2)
		{
			T_flow += (assign_lambda_eff(HostArraysPtr, local + def.locNx, def) * HostArraysPtr.T[local + def.locNx]
			- 2 * assign_lambda_eff(HostArraysPtr, local, def) * HostArraysPtr.T[local]
			+ assign_lambda_eff(HostArraysPtr, local - def.locNx, def) * HostArraysPtr.T[local - def.locNx]) / ((def.hy) * (def.hy));
		}

		if ((def.locNz) > 2)
		{
			T_flow = (assign_lambda_eff(HostArraysPtr, local + (def.locNx) * (def.locNy), def) * HostArraysPtr.T[local + (def.locNx) * (def.locNy)]
			- 2 * assign_lambda_eff(HostArraysPtr, local, def) * HostArraysPtr.T[local]
			+ assign_lambda_eff(HostArraysPtr, local - (def.locNx) * (def.locNy), def) * HostArraysPtr.T[local - (def.locNx) * (def.locNy)]) / ((def.hz) * (def.hz));
		}

		test_u(T_flow, __FILE__, __LINE__);
		return T_flow;
	}
	else
		return 0;
}

// Расчет потока энергии в точке
double assign_E_flow (ptr_Arrays HostArraysPtr, int i, int j, int k,  consts def)
{
	if ((i != 0) && (i != (def.locNx) - 1) && (j != 0) && (j != (def.locNy) - 1) && (((k != 0) && (k != (def.locNz) - 1)) || ((def.locNz) < 2)))
	{
		double E_flow = 0;
		int local=i + j * (def.locNx) + k * (def.locNx) * (def.locNy);

		if ((def.locNz) > 2)
		{
			E_flow += (HostArraysPtr.ro_w[local + 1] * HostArraysPtr.H_w[local + 1] * HostArraysPtr.ux_w[local + 1]
				- HostArraysPtr.ro_w[local - 1] * HostArraysPtr.H_w[local - 1] * HostArraysPtr.ux_w[local - 1]
				+ HostArraysPtr.ro_n[local + 1] * HostArraysPtr.H_n[local + 1] * HostArraysPtr.ux_n[local + 1]
				- HostArraysPtr.ro_n[local - 1] * HostArraysPtr.H_n[local - 1] * HostArraysPtr.ux_n[local - 1]
				+ HostArraysPtr.ro_g[local + 1] * HostArraysPtr.H_g[local + 1] * HostArraysPtr.ux_g[local + 1]
				- HostArraysPtr.ro_g[local - 1] * HostArraysPtr.H_g[local - 1] * HostArraysPtr.ux_g[local - 1]
				) / (2. * (def.hx));
		}
		if ((def.locNz) > 2)
		{
			E_flow += (HostArraysPtr.ro_w[local + def.locNx] * HostArraysPtr.H_w[local + def.locNx] * HostArraysPtr.uy_w[local + def.locNx]
			- HostArraysPtr.ro_w[local - def.locNx] * HostArraysPtr.H_w[local - def.locNx] * HostArraysPtr.uy_w[local - def.locNx]
			+ HostArraysPtr.ro_n[local + def.locNx] * HostArraysPtr.H_n[local + def.locNx] * HostArraysPtr.uy_n[local + def.locNx]
			- HostArraysPtr.ro_n[local - def.locNx] * HostArraysPtr.H_n[local - def.locNx] * HostArraysPtr.uy_n[local - def.locNx]
			+ HostArraysPtr.ro_g[local + def.locNx] * HostArraysPtr.H_g[local + def.locNx] * HostArraysPtr.uy_g[local + def.locNx]
			- HostArraysPtr.ro_g[local - def.locNx] * HostArraysPtr.H_g[local - def.locNx] * HostArraysPtr.uy_g[local - def.locNx]
			)/ (2. * (def.hy));	
		}

		if ((def.locNz) > 2)
		{
			E_flow += (HostArraysPtr.ro_w[local + (def.locNx) * (def.locNy)] * HostArraysPtr.H_w[local + (def.locNx) * (def.locNy)] * HostArraysPtr.uy_w[local + (def.locNx) * (def.locNy)]
			- HostArraysPtr.ro_w[local - (def.locNx) * (def.locNy)] * HostArraysPtr.H_w[local - (def.locNx) * (def.locNy)] * HostArraysPtr.uy_w[local - (def.locNx) * (def.locNy)]
			+ HostArraysPtr.ro_n[local + (def.locNx) * (def.locNy)] * HostArraysPtr.H_n[local + (def.locNx) * (def.locNy)] * HostArraysPtr.uy_n[local + (def.locNx) * (def.locNy)]
			- HostArraysPtr.ro_n[local - (def.locNx) * (def.locNy)] * HostArraysPtr.H_n[local - (def.locNx) * (def.locNy)] * HostArraysPtr.uy_n[local - (def.locNx) * (def.locNy)]
			+ HostArraysPtr.ro_g[local + (def.locNx) * (def.locNy)] * HostArraysPtr.H_g[local + (def.locNx) * (def.locNy)] * HostArraysPtr.uy_g[local + (def.locNx) * (def.locNy)]
			- HostArraysPtr.ro_g[local - (def.locNx) * (def.locNy)] * HostArraysPtr.H_g[local - (def.locNx) * (def.locNy)] * HostArraysPtr.uy_g[local - (def.locNx) * (def.locNy)]
			)/ (2. * (def.hz));	
		}

		test_u(E_flow, __FILE__, __LINE__);
		return E_flow;
	}
	else
		return 0;
}

// Расчет внутренней энергии всей системы в точке
void assign_E_current (ptr_Arrays HostArraysPtr, int local, consts def)
{
	HostArraysPtr.E[local] = (HostArraysPtr.m[local] * (HostArraysPtr.S_w[local] * (HostArraysPtr.ro_w[local] * HostArraysPtr.H_w[local] - HostArraysPtr.P_w[local])
		+ HostArraysPtr.S_n[local] * (HostArraysPtr.ro_n[local] * HostArraysPtr.H_n[local] - HostArraysPtr.P_n[local])
		+ HostArraysPtr.S_g[local] * (HostArraysPtr.ro_g[local] * HostArraysPtr.H_g[local] - HostArraysPtr.P_g[local])) 
		+ (1. - HostArraysPtr.m[local]) * (ro_r * HostArraysPtr.H_r[local] - HostArraysPtr.P_w[local]));

	test_nan(HostArraysPtr.E[local], __FILE__, __LINE__);
}

// Расчет внутренней энергии всей системы в точке на следующем шаге по времени
void assign_E_current (ptr_Arrays HostArraysPtr, int i, int j, int k, consts def)
{
	if ((i != 0) && (i != (def.locNx) - 1) && (j != 0) && (j != (def.locNy) - 1) && (((k != 0) && (k != (def.locNz) - 1)) || ((def.locNz) < 2)))
	{
		double Q_hw = 0, Q_hr = 0; // источниковые члены

		int local=i + j * (def.locNx) + k * (def.locNx) * (def.locNy);

		HostArraysPtr.E_new[local] = HostArraysPtr.E[local] + (def.dt) * (assign_T_flow(HostArraysPtr, i, j, k, def) + Q_hw + Q_hr - assign_E_flow(HostArraysPtr, i, j, k, def));

		test_nan(HostArraysPtr.E_new[local], __FILE__, __LINE__);
	}
}

// Задание граничных условий на температуру
void Border_T(ptr_Arrays HostArraysPtr, int i, int j, int k, consts def)
{
	if ((i == 0) || (i == (def.locNx) - 1) || (j == 0) || (j == (def.locNy) - 1) || (((k == 0) || (k == (def.locNz) - 1)) && ((def.locNz) >= 2)))
	{
		int local1 = set_boundary_basic_coordinate(i, j, k, def);
		int local = i + j * (def.locNx) + k * (def.locNx) * (def.locNy);

		// Будем считать границы области не теплопроводящими
		HostArraysPtr.T[local] = HostArraysPtr.T[local1];

		test_positive(HostArraysPtr.T[local], __FILE__, __LINE__);
	}
}