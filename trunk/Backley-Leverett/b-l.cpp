#include "../defines.h"

void assign_P2_Xi1_Xi2(ptr_Arrays HostArraysPtr, int i, int j, int k, int localNx, consts def)
{
	int media = HostArraysPtr.media[i+j*localNx+k*localNx*(def.Ny)];
	double S_e = (1. - HostArraysPtr.S2[i+j*localNx+k*localNx*(def.Ny)] - def.S_wr[media]) / (1. - def.S_wr[media]);
	double k1 = pow(S_e, (2. + 3. * def.lambda[media]) / def.lambda[media]);
	double k2 = (1. - S_e) * (1. - S_e) * (1 - pow(S_e, (2. + def.lambda[media]) / def.lambda[media]));
	double P_k = def.P_d[media] * pow((1. - HostArraysPtr.S2[i+j*localNx+k*localNx*(def.Ny)] - def.S_wr[media]) / (1. - def.S_wr[media]), -1. / def.lambda[media]);

	HostArraysPtr.P2[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.P1[i+j*localNx+k*localNx*(def.Ny)] + P_k;
	HostArraysPtr.Xi1[i+j*localNx+k*localNx*(def.Ny)] = -1. * def.K[media] * k1 / mu1;
	HostArraysPtr.Xi2[i+j*localNx+k*localNx*(def.Ny)] = -1. * def.K[media] * k2 / mu2;
}