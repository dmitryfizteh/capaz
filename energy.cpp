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
double �_w (double T, consts def)
{
	return c0_w - C_w * (T - T_0) + C_w2 * (T - T_0) * (T - T_0);
}

double �_n (double T, consts def)
{
	return c0_n + C_n * (T - T_0);
}

double �_g (double T, consts def)
{
	return c0_g + C_g * (T - T_0);
}

double �_r (double T, consts def)
{
	return c0_r + C_r * (T - T_0);
}

// ������������ ���������������� ��� water, napl, gas and rock
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

// ����������� ����������� ���������������� � ����� (����� �������������� ��� ������� ��������� ������)
double assign_lambda_eff (ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, int local, consts def)
{
	return HostArraysPtr.m[local] * (HostArraysPtr.S_w[local] * lambda_w (HostArraysPtr.T[local], def)
		+ HostArraysPtr.S_n[local] * lambda_n (HostArraysPtr.T[local], def)
		+ HostArraysPtr.S_g[local] * lambda_g (HostArraysPtr.T[local], def)) 
		+ (1. - HostArraysPtr.m[local]) * lambda_r (HostArraysPtr.T[local], def);
}

// ������ ��������� �� ����������� � ������������
// !!! ����������, ����� �������� ��� ������ �������� � ��������� ��������� �� �������, ����� �� ���������� ���� � �� ��
double assign_H_w (double T, consts def)
{
	double integral = 0, sum = 0, h_temp;
	int N_temp = 50;

	h_temp = (T - T_0) / N_temp;
	

	integral += �_w(T_0, def);
	integral += �_w(T, def);

	for(int i = 2; i < N_temp; i+=2)
		sum += �_w(T_0 + i * h_temp, def);

	sum *= 2;
	integral += sum;
	sum = 0;

	for(int i = 1; i < N_temp; i+=2)
		sum += �_w(T_0 + i * h_temp, def);

	sum *= 4;
	integral += sum;

	h_temp /= 3;
	integral *= h_temp;

	return integral;
}

double assign_H_n (double T, consts def)
{
	double integral = 0, sum = 0, h_temp;
	int N_temp = 50;

	h_temp = (T - T_0) / N_temp;


	integral += �_n(T_0, def);
	integral += �_n(T, def);

	for(int i = 2; i < N_temp; i+=2)
		sum += �_n(T_0 + i * h_temp, def);

	sum *= 2;
	integral += sum;
	sum = 0;

	for(int i = 1; i < N_temp; i+=2)
		sum += �_n(T_0 + i * h_temp, def);

	sum *= 4;
	integral += sum;

	h_temp /= 3;
	integral *= h_temp;

	return integral;
}

double assign_H_g (double T, consts def)
{
	double integral = 0, sum = 0, h_temp;
	int N_temp = 50;

	h_temp = (T - T_0) / N_temp;


	integral += �_g(T_0, def);
	integral += �_g(T, def);

	for(int i = 2; i < N_temp; i+=2)
		sum += �_g(T_0 + i * h_temp, def);

	sum *= 2;
	integral += sum;
	sum = 0;

	for(int i = 1; i < N_temp; i+=2)
		sum += �_g(T_0 + i * h_temp, def);

	sum *= 4;
	integral += sum;

	h_temp /= 3;
	integral *= h_temp;

	return integral;
}

void assign_H (ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, int local, consts def)
{
	assign_H_w (HostArraysPtr.T[local], def);
	assign_H_n (HostArraysPtr.T[local], def);
	assign_H_g (HostArraysPtr.T[local], def);
}

// ������������ �������� ��� water, napl, gas and rock

// ������ ������ ������� � �����

// ������ ��������� ������ � �����
