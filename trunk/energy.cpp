#include "defines.h"

// !!! ����� ����� ����� ������� � ��������� ��������
double T_0 = 273.; // �
double ro_r = 2000.; // ��/�^3
double lambda0_w = 0.553; // ��/(�*�)
double lambda0_n = 0.14;
double lambda0_g = 0.0237;
double lambda0_r = 1.;
double lambdaA_w = 3E-3; // 1/K
double lambdaA_n = 1E-3;
double lambdaA_g = 0.82;
double c0_w = 4.194E3; // ��/(��*�)
double c0_n = 1.7E3;
double c0_g = 1E3;
double c0_r = 0.8E3;
double C_w = 1.15;
double C_w2 = 0.015;
double C_n = 3.4;
double C_g = 0.119;
double C_r = 0.75;

// ������������ �������� ������������� ��� ���������� ��������  ��� water, napl, gas and rock � ��/(�*�)
double �_w (ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, int local, consts def)
{
	return c0_w - C_w * (HostArraysPtr.T[local] - T_0) + C_w2 * (HostArraysPtr.T[local] - T_0) * (HostArraysPtr.T[local] - T_0);
}

double �_n (ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, int local, consts def)
{
	return c0_n + C_n * (HostArraysPtr.T[local] - T_0);
}

double �_g (ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, int local, consts def)
{
	return c0_g + C_g * (HostArraysPtr.T[local] - T_0);
}

double �_r (ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, int local, consts def)
{
	return c0_r + C_r * (HostArraysPtr.T[local] - T_0);
}

// ������������ ���������������� ��� water, napl, gas and rock
double lambda_w (ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, int local, consts def)
{
	return lambda0_w * (1 - lambdaA_w * (HostArraysPtr.T[local] - T_0));
}

double lambda_n (ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, int local, consts def)
{
	return lambda0_n * (1 - lambdaA_n * (HostArraysPtr.T[local] - T_0));
}

double lambda_g (ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, int local, consts def)
{
	return lambda0_g * pow((HostArraysPtr.T[local] / T_0), lambdaA_g);
}

double lambda_r (ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, int local, consts def)
{
	return lambda0_r;
}

// ����������� ����������� ���������������� � ����� (����� �������������� ��� ������� ��������� ������)
double assign_lambda_eff (ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, int local, consts def)
{
	return HostArraysPtr.m[local] * (HostArraysPtr.S_w[local] * lambda_w (HostArraysPtr, DevArraysPtr, local, def)
		+ HostArraysPtr.S_n[local] * lambda_n (HostArraysPtr, DevArraysPtr, local, def)
		+ HostArraysPtr.S_g[local] * lambda_g (HostArraysPtr, DevArraysPtr, local, def)) 
		+ (1. - HostArraysPtr.m[local]) * lambda_r (HostArraysPtr, DevArraysPtr, local, def);
}

// ������ ��������� �� ����������� � ������������

// ������������ �������� ��� water, napl, gas and rock

// ������ ������ ������� � �����

// ������ ��������� ������ � �����
