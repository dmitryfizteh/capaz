#include "../defines.h"

void assign_Pn_Xi(ptr_Arrays HostArraysPtr, int i, int j, int k, int localNx, consts def)
{
	int media = HostArraysPtr.media[i+j*localNx+k*localNx*(def.Ny)];

	double S_w_e = (HostArraysPtr.S_w[i+j*localNx+k*localNx*(def.Ny)] - def.S_wr[media]) / (1. - def.S_wr[media]-def.S_nr[media]-def.S_gr[media]);
	double S_n_e = (HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)] - def.S_nr[media]) / (1. - def.S_wr[media]-def.S_nr[media]-def.S_gr[media]);
	double S_g_e = (HostArraysPtr.S_g[i+j*localNx+k*localNx*(def.Ny)] - def.S_gr[media]) / (1. - def.S_wr[media]-def.S_nr[media]-def.S_gr[media]);
	double k_w = pow(S_w_e,0.5)*pow(1.-pow(1.-pow(S_w_e,def.lambda[media]/(def.lambda[media]-1.)),(def.lambda[media]-1.)/def.lambda[media]),2.);
	double k_g = pow(S_g_e,0.5)*pow(1.-pow(1.-S_g_e,def.lambda[media]/(def.lambda[media]-1.)),2.*(def.lambda[media]-1.)/def.lambda[media]);
	double k_n_w = pow(1.-S_w_e,0.5)*pow(1.-pow(S_w_e,def.lambda[media]/(def.lambda[media]-1.)),2.*(def.lambda[media]-1.)/def.lambda[media]);
	double k_n_g = pow(S_n_e,0.5)*pow(1.-pow(1.-pow(S_n_e,def.lambda[media]/(def.lambda[media]-1.)),(def.lambda[media]-1.)/def.lambda[media]),2.);
	double k_n = S_n_e*k_n_w*k_n_g/(1-S_w_e)/(S_w_e+S_n_e);



	HostArraysPtr.P_n[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)];
	HostArraysPtr.Xi_w[i+j*localNx+k*localNx*(def.Ny)] = -1. * (def.K[media]) * k_w / def.mu_w;
	HostArraysPtr.Xi_n[i+j*localNx+k*localNx*(def.Ny)] = -1. * (def.K[media]) * k_n / def.mu_n;
}