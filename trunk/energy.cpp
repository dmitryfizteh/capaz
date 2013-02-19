#include "defines.h"

// !!! ����� ����� ����� ������� � ��������� ��������
double T_0 = 300; // �
double ro_r = 2000; // ��/�^3

// ������������ �������� ������������� ��� ���������� ��������  ��� water, napl, gas and rock � ��/(�*�)
double �_w (ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, int local, consts def)
{
	return 0;
}

double �_n (ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, int local, consts def)
{
	return 0;
}

double �_g (ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, int local, consts def)
{
	return 0;
}

double �_r (ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, int local, consts def)
{
	return 0.8 + 7.5E-4 * (HostArraysPtr.T[local] - 273.);
}

// ������������ ���������������� ��� water, napl, gas and rock
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
