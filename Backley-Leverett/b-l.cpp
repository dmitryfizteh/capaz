#include "../defines.h"
#include "b-l.h"

void assign_P_Xi(ptr_Arrays HostArraysPtr, int i, int j, int k, int localNx, consts def)
{
	int media = HostArraysPtr.media[i+j*localNx+k*localNx*(def.Ny)];
	double S_e = (1. - HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)] - def.S_wr[media]) / (1. - def.S_wr[media]);
	double k_w = pow(S_e, (2. + 3. * (def.lambda[media])) / def.lambda[media]);
	double k_n = (1. - S_e) * (1. - S_e) * (1 - pow(S_e, (2. + def.lambda[media]) / def.lambda[media]));
	
	HostArraysPtr.P_n[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)];
	HostArraysPtr.Xi_w[i+j*localNx+k*localNx*(def.Ny)] = -1. * (def.K[media]) * k_w / def.mu_w;
	HostArraysPtr.Xi_n[i+j*localNx+k*localNx*(def.Ny)] = -1. * (def.K[media]) * k_n / def.mu_n;

	test_nan(HostArraysPtr.P_n[i+j*localNx+k*localNx*(def.Ny)], __FILE__, __LINE__);
	test_nan(HostArraysPtr.Xi_w[i+j*localNx+k*localNx*(def.Ny)], __FILE__, __LINE__);
	test_nan(HostArraysPtr.Xi_n[i+j*localNx+k*localNx*(def.Ny)], __FILE__, __LINE__);
}

// Переписать!
void Newton(ptr_Arrays HostArraysPtr, int i, int j, int k, int localNx, consts def)
{
	if ((i!=0) && (i!=localNx-1) && (j!=0) && (j!=(def.Ny)-1) && (((k!=0) && (k!=(def.Nz)-1)) || ((def.Nz)<2)))
	{
		int media = HostArraysPtr.media[i+j*localNx+k*localNx*(def.Ny)];
		double S_e, P_k, AAA, F1, F2, PkS, F1P, F2P, F1S, F2S, det;

		for (int w=1;w<=def.newton_iterations;w++)
		{
			S_e = (1 - HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)] - def.S_wr[media]) / (1 - def.S_wr[media]);
			P_k=0;
			AAA = pow(S_e, (((-1.) / def.lambda[media]) - 1.));
			F1 = def.ro0_w * (1. + (def.beta_w) * (HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)] - def.P_atm)) * (1. - HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)]) - HostArraysPtr.roS_w[i+j*localNx+k*localNx*(def.Ny)];
			F2 = def.ro0_n * (1. + (def.beta_n) * (HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)] + P_k - def.P_atm)) * HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)] - HostArraysPtr.roS_n[i+j*localNx+k*localNx*(def.Ny)];
			PkS=0;
			F1P = def.ro0_w * (def.beta_w) * (1 - HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)]);
			F2P = def.ro0_n * (def.beta_n) * HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)];
			F1S = (-1) * (def.ro0_w) * (1 + (def.beta_w) * (HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)] - def.P_atm));
			F2S = def.ro0_n * (1 + (def.beta_n) * (HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)] + P_k - def.P_atm + (HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)] * PkS)));

			det = F1P * F2S - F1S * F2P;

			HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)] - (1 / det) * (F2S * F1 - F1S * F2);
			HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)] - (1 / det) * (F1P * F2 - F2P * F1);

			test_nan(HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)], __FILE__, __LINE__);
			test_nan(HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)], __FILE__, __LINE__);
		}  
	}
}

void Border_S(ptr_Arrays HostArraysPtr, int i, int j, int k, int localNx, int rank, int size, consts def)
{
	if ((i == 0) && ((def.Nx)>2))
	{
		HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.S_n[i+1+j*localNx+k*localNx*(def.Ny)];
		test_nan(HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)], __FILE__, __LINE__);
		return;
	}

	if ((i == localNx - 1) && ((def.Nx)>2))
	{
		HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.S_n[i-1+j*localNx+k*localNx*(def.Ny)];
		test_nan(HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)], __FILE__, __LINE__);
		return;
	}

	if ((j == (def.Ny) - 1) && ((def.Ny)>2))
	{
		HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.S_n[i+(j-1)*localNx+k*localNx*(def.Ny)];
		test_nan(HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)], __FILE__, __LINE__);
		return;
	}

	if ((j==0) && ((def.Ny)>2))
	{
		int I=i_to_I(i,rank,size, def);
		if ((I>=(def.Nx)/2-(def.source)) && (I<=(def.Nx)/2+(def.source)) && (k>=(def.Nz)/2-(def.source)) && (k<=(def.Nz)/2+(def.source)))
			HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)] = def.S_n_gr;
		else
			HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)] = 0;
		test_nan(HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)], __FILE__, __LINE__);
	}

	if ((k == 0) && ((def.Nz)>2))
	{
		HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.S_n[i+j*localNx+(k+1)*localNx*(def.Ny)];
		test_nan(HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)], __FILE__, __LINE__);
		return;
	}

	if ((k == (def.Nz) - 1) && ((def.Nz)>2))
	{
		HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.S_n[i+j*localNx+(k-1)*localNx*(def.Ny)];
		test_nan(HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)], __FILE__, __LINE__);
		return;
	}
}

void Border_P(ptr_Arrays HostArraysPtr, int i, int j, int k, int localNx, consts def)
{
	if ((i == 0) && ((def.Nx)>2))
	{
		HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.P_w[i+1+j*localNx+k*localNx*(def.Ny)]; 
		test_nan(HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)], __FILE__, __LINE__);
		return;
	}

	if ((i == localNx - 1) && ((def.Nx)>2))
	{
		HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.P_w[i-1+j*localNx+k*localNx*(def.Ny)];
		test_nan(HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)], __FILE__, __LINE__);
		return;
	}

	if ((j == (def.Ny) - 1) && ((def.Ny)>2))
	{
		HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.P_w[i+(j-1)*localNx+k*localNx*(def.Ny)] + HostArraysPtr.ro_w[i+localNx*1] * (def.g_const) * (def.hy);
		test_nan(HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)], __FILE__, __LINE__);
		return;
	}

	if ((j==0) && ((def.Ny)>2))
	{
		HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)] = def.P_atm;
		test_nan(HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)], __FILE__, __LINE__);
		return;
	}

	if ((k == 0) && ((def.Nz)>2))
	{
		HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.P_w[i+j*localNx+(k+1)*localNx*(def.Ny)]; 
		test_nan(HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)], __FILE__, __LINE__);
		return;
	}

	if ((k == (def.Nz) - 1) && ((def.Nz)>2))
	{
		HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.P_w[i+j*localNx+(k-1)*localNx*(def.Ny)];
		test_nan(HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)], __FILE__, __LINE__);
		return;
	}
}