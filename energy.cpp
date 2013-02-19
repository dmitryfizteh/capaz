#include "defines.h"

// !!! Нужно потом будет вынести в структуру констант
double T_0 = 300; // К
double ro_r = 2000; // кг/м^3

// Коэффициенты удельных теплоемкостей при постоянном давлении  для water, napl, gas and rock в Вт/(м*К)
double с_w (ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, int local, consts def)
{
	return 0;
}

double с_n (ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, int local, consts def)
{
	return 0;
}

double с_g (ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, int local, consts def)
{
	return 0;
}

double с_r (ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, int local, consts def)
{
	return 0.8 + 7.5E-4 * (HostArraysPtr.T[local] - 273.);
}

// Коэффициенты теплопроводности для water, napl, gas and rock
double lambda_w (ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, int local, consts def)
{
	return 0;
}

double lambda_n (ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, int local, consts def)
{
	return 0;
}

double lambda_g (ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, int local, consts def)
{
	return 0;
}

double lambda_r (ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, int local, consts def)
{
	return 0;
}

// Эффективный коэффициент теплопроводности в точке (будет использоваться при расчете теплового потока)
double assign_lambda_eff (ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, int local, consts def)
{
	return HostArraysPtr.m[local] * (HostArraysPtr.S_w[local] * lambda_w (HostArraysPtr, DevArraysPtr, local, def)
		+ HostArraysPtr.S_n[local] * lambda_n (HostArraysPtr, DevArraysPtr, local, def)
		+ HostArraysPtr.S_g[local] * lambda_g (HostArraysPtr, DevArraysPtr, local, def)) 
		+ (1. - HostArraysPtr.m[local]) * lambda_r (HostArraysPtr, DevArraysPtr, local, def);
}

// Расчет энтальпии по температуре и теплоемкости

// Коэффициенты вязкости для water, napl, gas and rock

// Расчет потока энергии в точке

// Расчет теплового потока в точке
