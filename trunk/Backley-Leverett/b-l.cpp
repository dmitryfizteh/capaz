#include "../defines.h"
#include "b-l.h"

void assign_P_Xi(ptr_Arrays HostArraysPtr, int i, int j, int k, localN locN, consts def)
{
	int media = HostArraysPtr.media[i+j*(locN.x)+k*(locN.x)*(locN.y)];
	double S_e = (1. - HostArraysPtr.S_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] - def.S_wr[media]) / (1. - def.S_wr[media]);
	//if (S_e<0)
	//	S_e=0;
	double k_w = pow(S_e, (2. + 3. * (def.lambda[media])) / def.lambda[media]);
	double k_n = (1. - S_e) * (1. - S_e) * (1 - pow(S_e, (2. + def.lambda[media]) / def.lambda[media]));

	HostArraysPtr.P_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] = HostArraysPtr.P_w[i+j*(locN.x)+k*(locN.x)*(locN.y)];
	HostArraysPtr.Xi_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] = -1. * (HostArraysPtr.K[i+j*(locN.x)+k*(locN.x)*(locN.y)]) * k_w / def.mu_w;
	HostArraysPtr.Xi_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] = -1. * (HostArraysPtr.K[i+j*(locN.x)+k*(locN.x)*(locN.y)]) * k_n / def.mu_n;

	test_S(S_e, __FILE__, __LINE__);
	test_positive(HostArraysPtr.P_n[i+j*(locN.x)+k*(locN.x)*(locN.y)], __FILE__, __LINE__);
	test_nan(HostArraysPtr.Xi_w[i+j*(locN.x)+k*(locN.x)*(locN.y)], __FILE__, __LINE__);
	test_nan(HostArraysPtr.Xi_n[i+j*(locN.x)+k*(locN.x)*(locN.y)], __FILE__, __LINE__);
}

// ѕереписать!
void Newton(ptr_Arrays HostArraysPtr, int i, int j, int k, localN locN, consts def)
{
	if ((i!=0) && (i!=(locN.x)-1) && (j!=0) && (j!=(locN.y)-1) && (((k!=0) && (k!=(locN.z)-1)) || ((locN.z)<2)))
	{
		int media = HostArraysPtr.media[i+j*(locN.x)+k*(locN.x)*(locN.y)];
		double S_e, F1, F2, F1P, F2P, F1S, F2S, det;

		for (int w=1;w<=def.newton_iterations;w++)
		{
			S_e = (1 - HostArraysPtr.S_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] - def.S_wr[media]) / (1 - def.S_wr[media]);

			F1 = def.ro0_w * (1. + (def.beta_w) * (HostArraysPtr.P_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] - def.P_atm)) * (1. - HostArraysPtr.S_n[i+j*(locN.x)+k*(locN.x)*(locN.y)]) - HostArraysPtr.roS_w[i+j*(locN.x)+k*(locN.x)*(locN.y)];
			F2 = def.ro0_n * (1. + (def.beta_n) * (HostArraysPtr.P_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] - def.P_atm)) * HostArraysPtr.S_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] - HostArraysPtr.roS_n[i+j*(locN.x)+k*(locN.x)*(locN.y)];
			
			F1P = def.ro0_w * (def.beta_w) * (1 - HostArraysPtr.S_n[i+j*(locN.x)+k*(locN.x)*(locN.y)]);
			F2P = def.ro0_n * (def.beta_n) * HostArraysPtr.S_n[i+j*(locN.x)+k*(locN.x)*(locN.y)];
			F1S = (-1) * (def.ro0_w) * (1 + (def.beta_w) * (HostArraysPtr.P_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] - def.P_atm));
			F2S = def.ro0_n * (1 + (def.beta_n) * (HostArraysPtr.P_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] - def.P_atm + (HostArraysPtr.S_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] )));

			det = F1P * F2S - F1S * F2P;

			HostArraysPtr.P_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] = HostArraysPtr.P_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] - (1 / det) * (F2S * F1 - F1S * F2);
			HostArraysPtr.S_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] = HostArraysPtr.S_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] - (1 / det) * (F1P * F2 - F2P * F1);

			test_positive(HostArraysPtr.P_w[i+j*(locN.x)+k*(locN.x)*(locN.y)], __FILE__, __LINE__);
			test_S(HostArraysPtr.S_n[i+j*(locN.x)+k*(locN.x)*(locN.y)], __FILE__, __LINE__);
		}  
	}
}

// «адание граничных условий с меньшим числом проверок, но с введением дополнительных переменных
void Border_S(ptr_Arrays HostArraysPtr, int i, int j, int k, localN locN, int rank, parts_sizes parts, consts def)
{
	int i1 = i, j1 = j, k1 = k;

	if(i == 0)
		i1 ++;
	if(i == (locN.x) - 1)
		i1 --;
	if(j == 0)
		j1 ++;
	if(j == (locN.y) - 1)
		j1 --;
	if((k == 0) && ((locN.z) > 2))
		k1 ++;
	if((k == (locN.z) - 1) && ((locN.z) > 2))
		k1 --;

	if((j != 0) || ((def.source) <= 0))
		HostArraysPtr.S_n[i + j * (locN.x) + k * (locN.x) * (locN.y)] = HostArraysPtr.S_n[i1 + j1 * (locN.x) + k1 * (locN.x) * (locN.y)];

	if((j == 0) && ((def.source) > 0))
	{
		int I=i_to_I(i, rank, parts, def);
		if ((I>=(def.Nx)/2-(def.source)) && (I<=(def.Nx)/2+(def.source)) && (k>=(def.Nz)/2-(def.source)) && (k<=(def.Nz)/2+(def.source)))
			HostArraysPtr.S_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] = def.S_n_gr;
		else
			//HostArraysPtr.S_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] = 0;
			HostArraysPtr.S_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] = HostArraysPtr.S_n[i+(j+1)*(locN.x)+k*(locN.x)*(locN.y)];
	}

	test_S(HostArraysPtr.S_n[i+j*(locN.x)+k*(locN.x)*(locN.y)], __FILE__, __LINE__);
}

void Border_P(ptr_Arrays HostArraysPtr, int i, int j, int k, localN locN, consts def)
{
	int i1 = i, j1 = j, k1 = k;

	if(i == 0)
		i1 ++;
	if(i == (locN.x) - 1)
		i1 --;
	if(j == 0)
		j1 ++;
	if(j == (locN.y) - 1)
		j1 --;
	if((k == 0) && ((locN.z) > 2))
		k1 ++;
	if((k == (locN.z) - 1) && ((locN.z) > 2))
		k1 --;

	if((j != 0) && (j != (locN.y) - 1))
		HostArraysPtr.P_w[i + j * (locN.x) + k * (locN.x) * (locN.y)] = HostArraysPtr.P_w[i1 + j1 * (locN.x) + k1 * (locN.x) * (locN.y)];
	else if(j == 0)
		HostArraysPtr.P_w[i + j * (locN.x) + k * (locN.x) * (locN.y)] = def.P_atm;
	else
		HostArraysPtr.P_w[i + j * (locN.x) + k * (locN.x) * (locN.y)] = HostArraysPtr.P_w[i1 + j1 * (locN.x) + k1 * (locN.x) * (locN.y)] + HostArraysPtr.ro_w[i1 + j1 * (locN.x) + k1 * (locN.x) * (locN.y)] * (def.g_const) * (def.hy);

	test_positive(HostArraysPtr.P_w[i+j*(locN.x)+k*(locN.x)*(locN.y)], __FILE__, __LINE__);
}

/*
void Border_S(ptr_Arrays HostArraysPtr, int i, int j, int k, int localNx, int rank, parts_sizes parts, consts def)
{
	if ((j == (def.Ny) - 1) && ((def.Ny)>2))
	{
		HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.S_n[i+(j-1)*localNx+k*localNx*(def.Ny)];
		test_S(HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)], __FILE__, __LINE__);
		return;
	}

	if ((j==0) && ((def.Ny)>2) && (j>0) && (j<def.Ny - 1))
	{
		int I=i_to_I(i,rank,parts, def);
		if ((I>=(def.Nx)/2-(def.source)) && (I<=(def.Nx)/2+(def.source)) && (k>=(def.Nz)/2-(def.source)) && (k<=(def.Nz)/2+(def.source)))
			HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)] = def.S_n_gr;
		else
			//HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)] = 0;
			HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.S_n[i+(j+1)*localNx+k*localNx*(def.Ny)];
		test_S(HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)], __FILE__, __LINE__);
		return;
	}

	if ((i == 0) && ((def.Nx)>2) && (j>0) && (j<def.Ny - 1))
	{
		HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.S_n[i+1+j*localNx+k*localNx*(def.Ny)];
		test_S(HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)], __FILE__, __LINE__);
		return;
	}

	if ((i == localNx - 1) && ((def.Nx)>2) && (j>0) && (j<def.Ny - 1))
	{
		HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.S_n[i-1+j*localNx+k*localNx*(def.Ny)];
		test_S(HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)], __FILE__, __LINE__);
		return;
	}

	if ((k == 0) && ((def.Nz)>2) && (j>0) && (j<def.Ny - 1))
	{
		HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.S_n[i+j*localNx+(k+1)*localNx*(def.Ny)];
		test_S(HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)], __FILE__, __LINE__);
		return;
	}

	if ((k == (def.Nz) - 1) && ((def.Nz)>2))
	{
		HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.S_n[i+j*localNx+(k-1)*localNx*(def.Ny)];
		test_S(HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)], __FILE__, __LINE__);
		return;
	}
}

void Border_P(ptr_Arrays HostArraysPtr, int i, int j, int k, int localNx, consts def)
{
	if ((j == (def.Ny) - 1) && ((def.Ny)>2))
	{
		HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.P_w[i+(j-1)*localNx+k*localNx*(def.Ny)] + HostArraysPtr.ro_w[i+(j-1)*localNx+k*localNx*(def.Ny)] * (def.g_const) * (def.hy); 
		test_positive(HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)], __FILE__, __LINE__);
		return;
	}

	if ((j==0) && ((def.Ny)>2))
	{
		HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)] = def.P_atm;
		test_positive(HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)], __FILE__, __LINE__);
		return;
	}

	if ((i == 0) && ((def.Nx)>2) && (j>0) && (j<def.Ny - 1))
	{
		HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.P_w[i+1+j*localNx+k*localNx*(def.Ny)]; 
		test_positive(HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)], __FILE__, __LINE__);
		return;
	}

	if ((i == localNx - 1) && ((def.Nx)>2) && (j>0) && (j<def.Ny - 1))
	{
		HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.P_w[i-1+j*localNx+k*localNx*(def.Ny)];
		test_positive(HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)], __FILE__, __LINE__);
		return;
	}

	if ((k == 0) && ((def.Nz)>2) && (j>0) && (j<def.Ny - 1))
	{
		HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.P_w[i+j*localNx+(k+1)*localNx*(def.Ny)]; 
		test_positive(HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)], __FILE__, __LINE__);
		return;
	}

	if ((k == (def.Nz) - 1) && ((def.Nz)>2) && (j>0) && (j<def.Ny - 1))
	{
		HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.P_w[i+j*localNx+(k-1)*localNx*(def.Ny)];
		test_positive(HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)], __FILE__, __LINE__);
		return;
	}
}
*/