#include "defines.h"

void ro_Pn_Xi_calculation(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, consts def, int localNx, int rank, int size, int blocksX, int blocksY, int blocksZ)
{
	for(int i=0;i<localNx;i++)
		for(int j=0;j<(def.Ny);j++)
			for(int k=0;k<(def.Nz);k++)
			{
				if(is_active_point(i, localNx, rank, size))
				{
					assign_ro(HostArraysPtr,i,j,k,localNx,def);
#ifdef THREE_PHASE
					assign_Pw_Pg_Xi(HostArraysPtr,i,j,k,localNx,def);
#else
					assign_Pn_Xi(HostArraysPtr,i,j,k,localNx,def);
#endif
				}
			}
}

void u_calculation(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, int localNx, int rank, int size, int blocksX, int blocksY, int blocksZ, consts def)
{
	for(int i=0;i<localNx;i++)
		for(int j=0;j<(def.Ny);j++)
			for(int k=0;k<(def.Nz);k++)
				if(is_active_point(i, localNx, rank, size))
					assign_u(HostArraysPtr,i,j,k,localNx,def);
}

void roS_calculation(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, consts def, double t, int localNx, int rank, int size, int blocksX, int blocksY, int blocksZ)
{
	for(int i=0;i<localNx;i++)
		for(int j=0;j<(def.Ny);j++)
			for(int k=0;k<(def.Nz);k++)
				if(is_active_point(i, localNx, rank, size))
				{
#ifdef NR
					assign_roS_nr(HostArraysPtr,t,i,j,k,localNx,def);
#else
					assign_roS(HostArraysPtr,t,i,j,k,localNx,def);
#endif
				}
}

void Pw_Sn_calculation(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, consts def, int localNx, int rank, int size, int blocksX, int blocksY, int blocksZ)
{
	for(int i=0;i<localNx;i++)
		for(int j=0;j<(def.Ny);j++)
			for(int k=0;k<(def.Nz);k++)
				if(is_active_point(i, localNx, rank, size))
					Newton(HostArraysPtr,i,j,k,localNx,def);
}


void boundary_conditions(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, int localNx, int rank, int size, int blocksX, int blocksY, int blocksZ, consts def)
{
	for(int i=0;i<localNx;i++)
		for(int j=0;j<(def.Ny);j++)
			for(int k=0;k<(def.Nz);k++)
				if(is_active_point(i, localNx, rank, size))
				{
#ifdef THREE_PHASE
					Border_Sw(HostArraysPtr,i,j,k,localNx,rank,size,def);
					Border_Sg(HostArraysPtr,i,j,k,localNx,rank,size,def);
					Border_Pn(HostArraysPtr,i,j,k,localNx,def);
#else
					Border_Sn(HostArraysPtr,i,j,k,localNx,rank,size,def);
					Border_Pw(HostArraysPtr,i,j,k,localNx,def);
#endif
				}
}	



void assign_ro(ptr_Arrays HostArraysPtr, int i, int j, int k, int localNx, consts def)
{
	HostArraysPtr.ro_w[i+j*localNx+k*localNx*(def.Ny)] = def.ro0_w * (1. + (def.beta_w) * (HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)] - def.P_atm));
	HostArraysPtr.ro_n[i+j*localNx+k*localNx*(def.Ny)] = def.ro0_n * (1. + (def.beta_n) * (HostArraysPtr.P_n[i+j*localNx+k*localNx*(def.Ny)] - def.P_atm));
#ifdef THREE_PHASE
	HostArraysPtr.ro_g[i+j*localNx+k*localNx*(def.Ny)] = def.ro0_g * (1. + (def.beta_g) * (HostArraysPtr.P_g[i+j*localNx+k*localNx*(def.Ny)] - def.P_atm));
#endif
}

// (!3D)
void assign_u(ptr_Arrays HostArraysPtr, int i, int j, int k, int localNx, consts def)
{
	//std::cout << i << " " << j << " " << k << "\n";
	if ((def.Nx)>2)
	{
		if (i == 0)
		{
			HostArraysPtr.ux_w[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.Xi_w[i+j*localNx+k*localNx*(def.Ny)] * (HostArraysPtr.P_w[i+1+j*localNx+k*localNx*(def.Ny)] - HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)]) / (def.hx);
			HostArraysPtr.ux_n[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.Xi_n[i+j*localNx+k*localNx*(def.Ny)] * (HostArraysPtr.P_n[i+1+j*localNx+k*localNx*(def.Ny)] - HostArraysPtr.P_n[i+j*localNx+k*localNx*(def.Ny)]) / (def.hx);
#ifdef THREE_PHASE
			HostArraysPtr.ux_g[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.Xi_g[i+j*localNx+k*localNx*(def.Ny)] * (HostArraysPtr.P_g[i+1+j*localNx+k*localNx*(def.Ny)] - HostArraysPtr.P_g[i+j*localNx+k*localNx*(def.Ny)]) / (def.hx);
#endif
		}
		else
		{
			if (i == localNx - 1)
			{
				HostArraysPtr.ux_w[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.Xi_w[i+j*localNx+k*localNx*(def.Ny)] * (HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)] - HostArraysPtr.P_w[i-1+j*localNx+k*localNx*(def.Ny)]) / (def.hx);
				HostArraysPtr.ux_n[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.Xi_n[i+j*localNx+k*localNx*(def.Ny)] * (HostArraysPtr.P_n[i+j*localNx+k*localNx*(def.Ny)] - HostArraysPtr.P_n[i-1+j*localNx+k*localNx*(def.Ny)]) / (def.hx);
#ifdef THREE_PHASE
				HostArraysPtr.ux_g[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.Xi_g[i+j*localNx+k*localNx*(def.Ny)] * (HostArraysPtr.P_g[i+j*localNx+k*localNx*(def.Ny)] - HostArraysPtr.P_g[i-1+j*localNx+k*localNx*(def.Ny)]) / (def.hx);
#endif
			}
			else
			{
				HostArraysPtr.ux_w[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.Xi_w[i+j*localNx+k*localNx*(def.Ny)] * ((HostArraysPtr.P_w[i+1+j*localNx+k*localNx*(def.Ny)] - HostArraysPtr.P_w[i-1+j*localNx+k*localNx*(def.Ny)]) / (2. * (def.hx)));
				HostArraysPtr.ux_n[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.Xi_n[i+j*localNx+k*localNx*(def.Ny)] * ((HostArraysPtr.P_n[i+1+j*localNx+k*localNx*(def.Ny)] - HostArraysPtr.P_n[i-1+j*localNx+k*localNx*(def.Ny)]) / (2. * (def.hx)));
#ifdef THREE_PHASE
				HostArraysPtr.ux_g[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.Xi_g[i+j*localNx+k*localNx*(def.Ny)] * ((HostArraysPtr.P_g[i+1+j*localNx+k*localNx*(def.Ny)] - HostArraysPtr.P_g[i-1+j*localNx+k*localNx*(def.Ny)]) / (2. * (def.hx)));
#endif			
			}
		}
	}
	else
	{
		HostArraysPtr.ux_w[i+j*localNx+k*localNx*(def.Ny)] = 0;
		HostArraysPtr.ux_n[i+j*localNx+k*localNx*(def.Ny)] = 0;
#ifdef THREE_PHASE
		HostArraysPtr.ux_g[i+j*localNx+k*localNx*(def.Ny)] = 0;
#endif
	}

	if ((def.Ny)>2)
	{
		if (j == 0)
		{
			HostArraysPtr.uy_w[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.Xi_w[i+j*localNx+k*localNx*(def.Ny)] * ((HostArraysPtr.P_w[i+(j+1)*localNx+k*localNx*(def.Ny)] - HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)]) / (def.hy) - HostArraysPtr.ro_w[i+j*localNx+k*localNx*(def.Ny)] * (def.g_const));
			HostArraysPtr.uy_n[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.Xi_n[i+j*localNx+k*localNx*(def.Ny)] * ((HostArraysPtr.P_n[i+(j+1)*localNx+k*localNx*(def.Ny)] - HostArraysPtr.P_n[i+j*localNx+k*localNx*(def.Ny)]) / (def.hy) - HostArraysPtr.ro_n[i+j*localNx+k*localNx*(def.Ny)] * (def.g_const));
#ifdef THREE_PHASE
			HostArraysPtr.uy_g[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.Xi_g[i+j*localNx+k*localNx*(def.Ny)] * ((HostArraysPtr.P_g[i+(j+1)*localNx+k*localNx*(def.Ny)] - HostArraysPtr.P_g[i+j*localNx+k*localNx*(def.Ny)]) / (def.hy) - HostArraysPtr.ro_g[i+j*localNx+k*localNx*(def.Ny)] * (def.g_const));
#endif
		}
		else
		{
			if (j == (def.Ny) - 1)
			{
				HostArraysPtr.uy_w[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.Xi_w[i+j*localNx+k*localNx*(def.Ny)] * ((HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)] - HostArraysPtr.P_w[i+(j-1)*localNx+k*localNx*(def.Ny)]) / (def.hy) - HostArraysPtr.ro_w[i+j*localNx+k*localNx*(def.Ny)] * (def.g_const));
				HostArraysPtr.uy_n[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.Xi_n[i+j*localNx+k*localNx*(def.Ny)] * ((HostArraysPtr.P_n[i+j*localNx+k*localNx*(def.Ny)] - HostArraysPtr.P_n[i+(j-1)*localNx+k*localNx*(def.Ny)]) / (def.hy) - HostArraysPtr.ro_n[i+j*localNx+k*localNx*(def.Ny)] * (def.g_const));
#ifdef THREE_PHASE
				HostArraysPtr.uy_g[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.Xi_g[i+j*localNx+k*localNx*(def.Ny)] * ((HostArraysPtr.P_g[i+j*localNx+k*localNx*(def.Ny)] - HostArraysPtr.P_g[i+(j-1)*localNx+k*localNx*(def.Ny)]) / (def.hy) - HostArraysPtr.ro_g[i+j*localNx+k*localNx*(def.Ny)] * (def.g_const));
#endif
			}
			else
			{
				HostArraysPtr.uy_w[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.Xi_w[i+j*localNx+k*localNx*(def.Ny)] * ((HostArraysPtr.P_w[i+(j+1)*localNx+k*localNx*(def.Ny)] - HostArraysPtr.P_w[i+(j-1)*localNx+k*localNx*(def.Ny)]) / (2. * (def.hy)) - HostArraysPtr.ro_w[i+j*localNx+k*localNx*(def.Ny)] * (def.g_const));
				HostArraysPtr.uy_n[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.Xi_n[i+j*localNx+k*localNx*(def.Ny)] * ((HostArraysPtr.P_n[i+(j+1)*localNx+k*localNx*(def.Ny)] - HostArraysPtr.P_n[i+(j-1)*localNx+k*localNx*(def.Ny)]) / (2. * (def.hy)) - HostArraysPtr.ro_n[i+j*localNx+k*localNx*(def.Ny)] * (def.g_const));
#ifdef THREE_PHASE
				HostArraysPtr.uy_g[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.Xi_g[i+j*localNx+k*localNx*(def.Ny)] * ((HostArraysPtr.P_g[i+(j+1)*localNx+k*localNx*(def.Ny)] - HostArraysPtr.P_g[i+(j-1)*localNx+k*localNx*(def.Ny)]) / (2. * (def.hy)) - HostArraysPtr.ro_g[i+j*localNx+k*localNx*(def.Ny)] * (def.g_const));
#endif
			}
		}
	}
	else
	{
		HostArraysPtr.uy_w[i+j*localNx+k*localNx*(def.Ny)] = 0;
		HostArraysPtr.uy_n[i+j*localNx+k*localNx*(def.Ny)] = 0;
#ifdef THREE_PHASE
		HostArraysPtr.uy_g[i+j*localNx+k*localNx*(def.Ny)] = 0;
#endif
	}

	if ((def.Nz)>2)
	{
		if (k == 0)
		{
			HostArraysPtr.uz_w[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.Xi_w[i+j*localNx+k*localNx*(def.Ny)] * (HostArraysPtr.P_w[i+localNx*j+(k+1)*localNx*(def.Ny)] - HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)]) / (def.hz);
			HostArraysPtr.uz_n[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.Xi_n[i+j*localNx+k*localNx*(def.Ny)] * (HostArraysPtr.P_n[i+localNx*j+(k+1)*localNx*(def.Ny)] - HostArraysPtr.P_n[i+j*localNx+k*localNx*(def.Ny)]) / (def.hz);
#ifdef THREE_PHASE
			HostArraysPtr.uz_g[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.Xi_g[i+j*localNx+k*localNx*(def.Ny)] * (HostArraysPtr.P_g[i+localNx*j+(k+1)*localNx*(def.Ny)] - HostArraysPtr.P_g[i+j*localNx+k*localNx*(def.Ny)]) / (def.hz);
#endif
		}
		else
		{
			if (k == (def.Nz) - 1)
			{
				HostArraysPtr.uz_w[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.Xi_w[i+j*localNx+k*localNx*(def.Ny)] * (HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)] - HostArraysPtr.P_w[i+localNx*j+(k-1)*localNx*(def.Ny)]) / (def.hz);
				HostArraysPtr.uz_n[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.Xi_n[i+j*localNx+k*localNx*(def.Ny)] * (HostArraysPtr.P_n[i+j*localNx+k*localNx*(def.Ny)] - HostArraysPtr.P_n[i+localNx*j+(k-1)*localNx*(def.Ny)]) / (def.hz);
#ifdef THREE_PHASE
				HostArraysPtr.uz_g[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.Xi_g[i+j*localNx+k*localNx*(def.Ny)] * (HostArraysPtr.P_g[i+j*localNx+k*localNx*(def.Ny)] - HostArraysPtr.P_g[i+localNx*j+(k-1)*localNx*(def.Ny)]) / (def.hz);
#endif
			}
			else
			{
				HostArraysPtr.uz_w[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.Xi_w[i+j*localNx+k*localNx*(def.Ny)] * ((HostArraysPtr.P_w[i+localNx*j+(k+1)*localNx*(def.Ny)] - HostArraysPtr.P_w[i+localNx*j+(k-1)*localNx*(def.Ny)]) / (2. * (def.hz)));
				HostArraysPtr.uz_n[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.Xi_n[i+j*localNx+k*localNx*(def.Ny)] * ((HostArraysPtr.P_n[i+localNx*j+(k+1)*localNx*(def.Ny)] - HostArraysPtr.P_n[i+localNx*j+(k-1)*localNx*(def.Ny)]) / (2. * (def.hz)));
#ifdef THREE_PHASE
				HostArraysPtr.uz_g[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.Xi_g[i+j*localNx+k*localNx*(def.Ny)] * ((HostArraysPtr.P_g[i+localNx*j+(k+1)*localNx*(def.Ny)] - HostArraysPtr.P_g[i+localNx*j+(k-1)*localNx*(def.Ny)]) / (2. * (def.hz)));
#endif
			}
		}
	}
	else
	{
		HostArraysPtr.uz_w[i+j*localNx+k*localNx*(def.Ny)] = 0;
		HostArraysPtr.uz_n[i+j*localNx+k*localNx*(def.Ny)] = 0;
#ifdef THREE_PHASE
		HostArraysPtr.uz_g[i+j*localNx+k*localNx*(def.Ny)] = 0;
#endif
	}
}

#ifdef THREE_PHASE
void assign_roS(ptr_Arrays HostArraysPtr, double t, int i, int j, int k, int localNx, consts def)
{
	if ((i!=0) && (i!=localNx-1) && (j!=0) && (j!=(def.Ny)-1) && (((k!=0) && (k!=(def.Nz)-1)) || ((def.Nz)<2)))
	{
		int media = HostArraysPtr.media[i+j*localNx+k*localNx*(def.Ny)];

		HostArraysPtr.roS_w[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.ro_w[i+j*localNx+k*localNx*(def.Ny)] * HostArraysPtr.S_w[i+j*localNx+k*localNx*(def.Ny)];
		HostArraysPtr.roS_n[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.ro_n[i+j*localNx+k*localNx*(def.Ny)] * HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)];
		HostArraysPtr.roS_g[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.ro_g[i+j*localNx+k*localNx*(def.Ny)] * HostArraysPtr.S_g[i+j*localNx+k*localNx*(def.Ny)];

		double divgrad1, divgrad2, divgrad3, Tx1, Ty1, Tz1, Tx2, Ty2, Tz2, Tx3, Ty3, Tz3, A1=0, A2=0, A3=0;

		if ((def.Nz)<2)
		{
			divgrad1=0;
			divgrad2=0;
			divgrad3=0;
			Tz1=0;
			Tz2=0;
			Tz3=0;
		}
		else
		{
			divgrad1 = (def.m[media] * (def.l_w) * (def.c) / 2.) * (HostArraysPtr.ro_w[i+j*localNx+(k+1)*localNx*(def.Ny)] * HostArraysPtr.S_w[i+j*localNx+(k+1)*localNx*(def.Ny)] - 2 * HostArraysPtr.ro_w[i+j*localNx+k*localNx*(def.Ny)] * HostArraysPtr.S_w[i+j*localNx+k*localNx*(def.Ny)] + HostArraysPtr.ro_w[i+j*localNx+(k-1)*localNx*(def.Ny)] * (HostArraysPtr.S_w[i+j*localNx+(k-1)*localNx*(def.Ny)])) / ((def.hz) * (def.hz));
			divgrad2 = (def.m[media] * (def.l_n) * (def.c) / 2.) * (HostArraysPtr.ro_n[i+j*localNx+(k+1)*localNx*(def.Ny)] * HostArraysPtr.S_n[i+j*localNx+(k+1)*localNx*(def.Ny)] - 2 * HostArraysPtr.ro_n[i+j*localNx+k*localNx*(def.Ny)] * HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)] + HostArraysPtr.ro_n[i+j*localNx+(k-1)*localNx*(def.Ny)] * (HostArraysPtr.S_n[i+j*localNx+(k-1)*localNx*(def.Ny)])) / ((def.hz) * (def.hz));
			
			Tz1 = (HostArraysPtr.ro_w[i+1+j*localNx+(k+1)*localNx*(def.Ny)] * HostArraysPtr.ux_w[i+1+j*localNx+(k+1)*localNx*(def.Ny)] - HostArraysPtr.ro_w[i+j*localNx+(k-1)*localNx*(def.Ny)] * HostArraysPtr.ux_w[i+j*localNx+(k-1)*localNx*(def.Ny)]) / (2. * (def.hz));
			Tz2 = (HostArraysPtr.ro_n[i+j*localNx+(k+1)*localNx*(def.Ny)] * HostArraysPtr.uy_n[i+j*localNx+(k+1)*localNx*(def.Ny)] - HostArraysPtr.ro_n[i+j*localNx+(k-1)*localNx*(def.Ny)] * HostArraysPtr.uy_n[i+j*localNx+(k-1)*localNx*(def.Ny)]) / (2. * (def.hz));
		}

		divgrad1 += (def.m[media] * (def.l_w) * (def.c) / 2.) *
			((HostArraysPtr.ro_w[i+1+j*localNx+k*localNx*(def.Ny)] * (1. - HostArraysPtr.S_n[i+1+j*localNx+k*localNx*(def.Ny)]) - 2 * HostArraysPtr.ro_w[i+j*localNx+k*localNx*(def.Ny)] * (1. - HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)]) + HostArraysPtr.ro_w[i-1+j*localNx+k*localNx*(def.Ny)] * (1. - HostArraysPtr.S_n[i-1+j*localNx+k*localNx*(def.Ny)])) / ((def.hx) * (def.hx)) +
			(HostArraysPtr.ro_w[i+(j+1)*localNx+k*localNx*(def.Ny)] * (1. - HostArraysPtr.S_n[i+(j+1)*localNx+k*localNx*(def.Ny)]) - 2 * HostArraysPtr.ro_w[i+j*localNx+k*localNx*(def.Ny)] * (1. - HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)]) + HostArraysPtr.ro_w[i+(j-1)*localNx+k*localNx*(def.Ny)] * (1. - HostArraysPtr.S_n[i+(j-1)*localNx+k*localNx*(def.Ny)])) / ((def.hy) * (def.hy)));

		divgrad2 += (def.m[media] * (def.l_n) * (def.c) / 2.) *
			((HostArraysPtr.ro_n[i+1+j*localNx+k*localNx*(def.Ny)] * HostArraysPtr.S_n[i+1+j*localNx+k*localNx*(def.Ny)] - 2 * HostArraysPtr.ro_n[i+j*localNx+k*localNx*(def.Ny)] * HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)] + HostArraysPtr.ro_n[i-1+j*localNx+k*localNx*(def.Ny)] * (HostArraysPtr.S_n[i-1+j*localNx+k*localNx*(def.Ny)])) / ((def.hx) * (def.hx)) +
			(HostArraysPtr.ro_n[i+(j+1)*localNx+k*localNx*(def.Ny)] * HostArraysPtr.S_n[i+(j+1)*localNx+k*localNx*(def.Ny)] - 2 * HostArraysPtr.ro_n[i+j*localNx+k*localNx*(def.Ny)] * HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)] + HostArraysPtr.ro_n[i+(j-1)*localNx+k*localNx*(def.Ny)] * (HostArraysPtr.S_n[i+(j-1)*localNx+k*localNx*(def.Ny)])) / ((def.hy) * (def.hy)));

		Tx1 = (HostArraysPtr.ro_w[i+1+j*localNx+k*localNx*(def.Ny)] * HostArraysPtr.ux_w[i+1+j*localNx+k*localNx*(def.Ny)] - HostArraysPtr.ro_w[i-1+j*localNx+k*localNx*(def.Ny)] * HostArraysPtr.ux_w[i-1+j*localNx+k*localNx*(def.Ny)]) / (2. * (def.hx));
		Ty1 = (HostArraysPtr.ro_w[i+(j+1)*localNx+k*localNx*(def.Ny)] * HostArraysPtr.uy_w[i+(j+1)*localNx+k*localNx*(def.Ny)] - HostArraysPtr.ro_w[i+(j-1)*localNx+k*localNx*(def.Ny)] * HostArraysPtr.uy_w[i+(j-1)*localNx+k*localNx*(def.Ny)]) / (2. * (def.hy));
		Tx2 = (HostArraysPtr.ro_n[i+1+j*localNx+k*localNx*(def.Ny)] * HostArraysPtr.ux_n[i+1+j*localNx+k*localNx*(def.Ny)] - HostArraysPtr.ro_n[i-1+j*localNx+k*localNx*(def.Ny)] * HostArraysPtr.ux_n[i-1+j*localNx+k*localNx*(def.Ny)]) / (2. * (def.hx));
		Ty2 = (HostArraysPtr.ro_n[i+(j+1)*localNx+k*localNx*(def.Ny)] * HostArraysPtr.uy_n[i+(j+1)*localNx+k*localNx*(def.Ny)] - HostArraysPtr.ro_n[i+(j-1)*localNx+k*localNx*(def.Ny)] * HostArraysPtr.uy_n[i+(j-1)*localNx+k*localNx*(def.Ny)]) / (2. * (def.hy));

		if (t < 2 * (def.dt))
		{
			A1 = HostArraysPtr.roS_w[i+j*localNx+k*localNx*(def.Ny)] + ((def.dt) / def.m[media]) * (divgrad1 - Tx1 - Ty1 - Tz1);
			A2 = HostArraysPtr.roS_n[i+j*localNx+k*localNx*(def.Ny)] + ((def.dt) / def.m[media]) * (divgrad2 - Tx2 - Ty2 - Tz2);
		}
		else
		{
			A1 = (2. * (def.dt) * (def.dt)) / (def.m[media] * ((def.dt) + 2. * (def.tau))) * (divgrad1 - Tx1 - Ty1 - Tz1 + (2. * HostArraysPtr.roS_w[i+j*localNx+k*localNx*(def.Ny)] * (def.m[media]) * (def.tau)) / ((def.dt) * (def.dt)) + HostArraysPtr.roS_w_old[i+j*localNx+k*localNx*(def.Ny)] * (def.m[media]) * ((def.dt) - 2. * (def.tau)) / (2. * (def.dt) * (def.dt)));
			A2 = (2. * (def.dt) * (def.dt)) / (def.m[media] * ((def.dt) + 2. * (def.tau))) * (divgrad2 - Tx2 - Ty2 - Tz2 + (2. * HostArraysPtr.roS_n[i+j*localNx+k*localNx*(def.Ny)] * (def.m[media]) * (def.tau)) / ((def.dt) * (def.dt)) + HostArraysPtr.roS_n_old[i+j*localNx+k*localNx*(def.Ny)] * (def.m[media]) * ((def.dt) - 2. * (def.tau)) / (2. * (def.dt) * (def.dt)));
		}
		HostArraysPtr.roS_w_old[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.roS_w[i+j*localNx+k*localNx*(def.Ny)];
		HostArraysPtr.roS_n_old[i+j*localNx+k*localNx*(def.Ny)]= HostArraysPtr.roS_n[i+j*localNx+k*localNx*(def.Ny)];
		HostArraysPtr.roS_w[i+j*localNx+k*localNx*(def.Ny)] = A1;
		HostArraysPtr.roS_n[i+j*localNx+k*localNx*(def.Ny)] = A2;
	}
}
#else
void assign_roS(ptr_Arrays HostArraysPtr, double t, int i, int j, int k, int localNx, consts def)
{
	if ((i!=0) && (i!=localNx-1) && (j!=0) && (j!=(def.Ny)-1) && (((k!=0) && (k!=(def.Nz)-1)) || ((def.Nz)<2)))
	{
		int media = HostArraysPtr.media[i+j*localNx+k*localNx*(def.Ny)];

		HostArraysPtr.roS_w[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.ro_w[i+j*localNx+k*localNx*(def.Ny)] * (1 - HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)]);
		HostArraysPtr.roS_n[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.ro_n[i+j*localNx+k*localNx*(def.Ny)] * HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)];

		double divgrad1, divgrad2, Tx1, Ty1, Tz1, Tx2, Ty2, Tz2, A1=0, A2=0;

		if ((def.Nz)<2)
		{
			divgrad1=0;
			divgrad2=0;
			Tz1=0;
			Tz2=0;
		}
		else
		{
			divgrad1 = (def.m[media] * (def.l_w) * (def.c) / 2.) * (HostArraysPtr.ro_w[i+j*localNx+(k+1)*localNx*(def.Ny)] * (1. - HostArraysPtr.S_n[i+j*localNx+(k+1)*localNx*(def.Ny)]) - 2 * HostArraysPtr.ro_w[i+j*localNx+k*localNx*(def.Ny)] * (1. - HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)]) + HostArraysPtr.ro_w[i+j*localNx+(k-1)*localNx*(def.Ny)] * (1. - HostArraysPtr.S_n[i+j*localNx+(k-1)*localNx*(def.Ny)])) / ((def.hz) * (def.hz));
			divgrad2 = (def.m[media] * (def.l_n) * (def.c) / 2.) * (HostArraysPtr.ro_n[i+j*localNx+(k+1)*localNx*(def.Ny)] * HostArraysPtr.S_n[i+j*localNx+(k+1)*localNx*(def.Ny)] - 2 * HostArraysPtr.ro_n[i+j*localNx+k*localNx*(def.Ny)] * HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)] + HostArraysPtr.ro_n[i+j*localNx+(k-1)*localNx*(def.Ny)] * (HostArraysPtr.S_n[i+j*localNx+(k-1)*localNx*(def.Ny)])) / ((def.hz) * (def.hz));
			Tz1 = (HostArraysPtr.ro_w[i+1+j*localNx+(k+1)*localNx*(def.Ny)] * HostArraysPtr.ux_w[i+1+j*localNx+(k+1)*localNx*(def.Ny)] - HostArraysPtr.ro_w[i+j*localNx+(k-1)*localNx*(def.Ny)] * HostArraysPtr.ux_w[i+j*localNx+(k-1)*localNx*(def.Ny)]) / (2. * (def.hz));
			Tz2 = (HostArraysPtr.ro_n[i+j*localNx+(k+1)*localNx*(def.Ny)] * HostArraysPtr.uy_n[i+j*localNx+(k+1)*localNx*(def.Ny)] - HostArraysPtr.ro_n[i+j*localNx+(k-1)*localNx*(def.Ny)] * HostArraysPtr.uy_n[i+j*localNx+(k-1)*localNx*(def.Ny)]) / (2. * (def.hz));
		}

		divgrad1 += (def.m[media] * (def.l_w) * (def.c) / 2.) *
			((HostArraysPtr.ro_w[i+1+j*localNx+k*localNx*(def.Ny)] * (1. - HostArraysPtr.S_n[i+1+j*localNx+k*localNx*(def.Ny)]) - 2 * HostArraysPtr.ro_w[i+j*localNx+k*localNx*(def.Ny)] * (1. - HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)]) + HostArraysPtr.ro_w[i-1+j*localNx+k*localNx*(def.Ny)] * (1. - HostArraysPtr.S_n[i-1+j*localNx+k*localNx*(def.Ny)])) / ((def.hx) * (def.hx)) +
			(HostArraysPtr.ro_w[i+(j+1)*localNx+k*localNx*(def.Ny)] * (1. - HostArraysPtr.S_n[i+(j+1)*localNx+k*localNx*(def.Ny)]) - 2 * HostArraysPtr.ro_w[i+j*localNx+k*localNx*(def.Ny)] * (1. - HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)]) + HostArraysPtr.ro_w[i+(j-1)*localNx+k*localNx*(def.Ny)] * (1. - HostArraysPtr.S_n[i+(j-1)*localNx+k*localNx*(def.Ny)])) / ((def.hy) * (def.hy)));

		divgrad2 += (def.m[media] * (def.l_n) * (def.c) / 2.) *
			((HostArraysPtr.ro_n[i+1+j*localNx+k*localNx*(def.Ny)] * HostArraysPtr.S_n[i+1+j*localNx+k*localNx*(def.Ny)] - 2 * HostArraysPtr.ro_n[i+j*localNx+k*localNx*(def.Ny)] * HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)] + HostArraysPtr.ro_n[i-1+j*localNx+k*localNx*(def.Ny)] * (HostArraysPtr.S_n[i-1+j*localNx+k*localNx*(def.Ny)])) / ((def.hx) * (def.hx)) +
			(HostArraysPtr.ro_n[i+(j+1)*localNx+k*localNx*(def.Ny)] * HostArraysPtr.S_n[i+(j+1)*localNx+k*localNx*(def.Ny)] - 2 * HostArraysPtr.ro_n[i+j*localNx+k*localNx*(def.Ny)] * HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)] + HostArraysPtr.ro_n[i+(j-1)*localNx+k*localNx*(def.Ny)] * (HostArraysPtr.S_n[i+(j-1)*localNx+k*localNx*(def.Ny)])) / ((def.hy) * (def.hy)));

		Tx1 = (HostArraysPtr.ro_w[i+1+j*localNx+k*localNx*(def.Ny)] * HostArraysPtr.ux_w[i+1+j*localNx+k*localNx*(def.Ny)] - HostArraysPtr.ro_w[i-1+j*localNx+k*localNx*(def.Ny)] * HostArraysPtr.ux_w[i-1+j*localNx+k*localNx*(def.Ny)]) / (2. * (def.hx));
		Ty1 = (HostArraysPtr.ro_w[i+(j+1)*localNx+k*localNx*(def.Ny)] * HostArraysPtr.uy_w[i+(j+1)*localNx+k*localNx*(def.Ny)] - HostArraysPtr.ro_w[i+(j-1)*localNx+k*localNx*(def.Ny)] * HostArraysPtr.uy_w[i+(j-1)*localNx+k*localNx*(def.Ny)]) / (2. * (def.hy));
		Tx2 = (HostArraysPtr.ro_n[i+1+j*localNx+k*localNx*(def.Ny)] * HostArraysPtr.ux_n[i+1+j*localNx+k*localNx*(def.Ny)] - HostArraysPtr.ro_n[i-1+j*localNx+k*localNx*(def.Ny)] * HostArraysPtr.ux_n[i-1+j*localNx+k*localNx*(def.Ny)]) / (2. * (def.hx));
		Ty2 = (HostArraysPtr.ro_n[i+(j+1)*localNx+k*localNx*(def.Ny)] * HostArraysPtr.uy_n[i+(j+1)*localNx+k*localNx*(def.Ny)] - HostArraysPtr.ro_n[i+(j-1)*localNx+k*localNx*(def.Ny)] * HostArraysPtr.uy_n[i+(j-1)*localNx+k*localNx*(def.Ny)]) / (2. * (def.hy));

		if (t < 2 * (def.dt))
		{
			A1 = HostArraysPtr.roS_w[i+j*localNx+k*localNx*(def.Ny)] + ((def.dt) / def.m[media]) * (divgrad1 - Tx1 - Ty1 - Tz1);
			A2 = HostArraysPtr.roS_n[i+j*localNx+k*localNx*(def.Ny)] + ((def.dt) / def.m[media]) * (divgrad2 - Tx2 - Ty2 - Tz2);
		}
		else
		{
			A1 = (2. * (def.dt) * (def.dt)) / (def.m[media] * ((def.dt) + 2. * (def.tau))) * (divgrad1 - Tx1 - Ty1 - Tz1 + (2. * HostArraysPtr.roS_w[i+j*localNx+k*localNx*(def.Ny)] * (def.m[media]) * (def.tau)) / ((def.dt) * (def.dt)) + HostArraysPtr.roS_w_old[i+j*localNx+k*localNx*(def.Ny)] * (def.m[media]) * ((def.dt) - 2. * (def.tau)) / (2. * (def.dt) * (def.dt)));
			A2 = (2. * (def.dt) * (def.dt)) / (def.m[media] * ((def.dt) + 2. * (def.tau))) * (divgrad2 - Tx2 - Ty2 - Tz2 + (2. * HostArraysPtr.roS_n[i+j*localNx+k*localNx*(def.Ny)] * (def.m[media]) * (def.tau)) / ((def.dt) * (def.dt)) + HostArraysPtr.roS_n_old[i+j*localNx+k*localNx*(def.Ny)] * (def.m[media]) * ((def.dt) - 2. * (def.tau)) / (2. * (def.dt) * (def.dt)));
		}
		HostArraysPtr.roS_w_old[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.roS_w[i+j*localNx+k*localNx*(def.Ny)];
		HostArraysPtr.roS_n_old[i+j*localNx+k*localNx*(def.Ny)]= HostArraysPtr.roS_n[i+j*localNx+k*localNx*(def.Ny)];
		HostArraysPtr.roS_w[i+j*localNx+k*localNx*(def.Ny)] = A1;
		HostArraysPtr.roS_n[i+j*localNx+k*localNx*(def.Ny)] = A2;
	}
}
#endif

void assign_roS_nr(ptr_Arrays HostArraysPtr, double t, int i, int j, int k, int localNx, consts def)
{
	if ((i!=0) && (i!=localNx-1) && (j!=0) && (j!=(def.Ny)-1) && (((k!=0) && (k!=(def.Nz)-1)) || ((def.Nz)<2)))
	{
		int media = HostArraysPtr.media[i+j*localNx+k*localNx*(def.Ny)];

		HostArraysPtr.roS_w[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.ro_w[i+j*localNx+k*localNx*(def.Ny)] * (1 - HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)]);
		HostArraysPtr.roS_n[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.ro_n[i+j*localNx+k*localNx*(def.Ny)] * HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)];

		double P1 = HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)];
		double P2 = HostArraysPtr.P_n[i+j*localNx+k*localNx*(def.Ny)];

		double x1, x2, y1, y2, z1, z2, f1, f2, f3, g1, g2, g3, A1=0, A2=0;

		if ((def.Nz)<2)
		{
			f3=0;
			g3=0;
		}
		else
		{
			z2 = -(HostArraysPtr.P_w[i+j*localNx+(k+1)*localNx*(def.Ny)] - P1)/def.hz;
			z1 = -(P1 - HostArraysPtr.P_w[i+j*localNx+(k-1)*localNx*(def.Ny)])/def.hz;

			f3 = (((z2 + abs(z2))/2.0 - (z1 - abs(z1))/2.0)*(-1) * HostArraysPtr.Xi_w[i+j*localNx+k*localNx*(def.Ny)] * HostArraysPtr.ro_w[i+j*localNx+k*localNx*(def.Ny)] -
                      (z1 + abs(z1))/2.0*(-1)* HostArraysPtr.Xi_w[i+j*localNx+(k-1)*localNx*(def.Ny)] * HostArraysPtr.ro_w[i+j*localNx+(k-1)*localNx*(def.Ny)] +
                      (z2 - abs(z2))/2.0*(-1)* HostArraysPtr.Xi_w[i+j*localNx+(k+1)*localNx*(def.Ny)] * HostArraysPtr.ro_w[i+j*localNx+(k+1)*localNx*(def.Ny)])/def.hz;

			z2 = -(HostArraysPtr.P_n[i+j*localNx+(k+1)*localNx*(def.Ny)] - P2)/def.hz;
			z1 = -(P2 - HostArraysPtr.P_n[i+j*localNx+(k-1)*localNx*(def.Ny)])/def.hz;

			g3 = (((z2 + abs(z2))/2.0 - (z1 - abs(z1))/2.0)*(-1) * HostArraysPtr.Xi_n[i+j*localNx+k*localNx*(def.Ny)] * HostArraysPtr.ro_n[i+j*localNx+k*localNx*(def.Ny)] -
                      (z1 + abs(z1))/2.0*(-1)* HostArraysPtr.Xi_n[i+j*localNx+(k-1)*localNx*(def.Ny)] * HostArraysPtr.ro_n[i+j*localNx+(k-1)*localNx*(def.Ny)] +
                      (z2 - abs(z2))/2.0*(-1)* HostArraysPtr.Xi_n[i+j*localNx+(k+1)*localNx*(def.Ny)] * HostArraysPtr.ro_n[i+j*localNx+(k+1)*localNx*(def.Ny)])/def.hz;
		}

		x2 = -(HostArraysPtr.P_w[i+1+j*localNx+k*localNx*(def.Ny)] - P1)/def.hx;
        x1 = -(P1 - HostArraysPtr.P_w[i-1+j*localNx+k*localNx*(def.Ny)])/def.hx;

        y2 = -(HostArraysPtr.P_w[i+(j+1)*localNx+k*localNx*(def.Ny)] - P1)/def.hy + HostArraysPtr.ro_w[i+j*localNx+k*localNx*(def.Ny)] * (def.g_const);
        y1 = -(P1 - HostArraysPtr.P_w[i+(j-1)*localNx+k*localNx*(def.Ny)])/def.hy + HostArraysPtr.ro_w[i+j*localNx+k*localNx*(def.Ny)] * (def.g_const);

        f1 = (((x2 + abs(x2))/2.0 - (x1 - abs(x1))/2.0)*(-1) * HostArraysPtr.Xi_w[i+j*localNx+k*localNx*(def.Ny)] * HostArraysPtr.ro_w[i+j*localNx+k*localNx*(def.Ny)] -
                (x1 + abs(x1))/2.0*(-1)* HostArraysPtr.Xi_w[i-1+j*localNx+k*localNx*(def.Ny)] * HostArraysPtr.ro_w[i-1+j*localNx+k*localNx*(def.Ny)] +
                (x2 - abs(x2))/2.0*(-1)* HostArraysPtr.Xi_w[i+1+j*localNx+k*localNx*(def.Ny)] * HostArraysPtr.ro_w[i+1+j*localNx+k*localNx*(def.Ny)])/def.hx;

        f2 = (((y2 + abs(y2))/2.0 - (y1 - abs(y1))/2.0)*(-1)* HostArraysPtr.Xi_w[i+j*localNx+k*localNx*(def.Ny)] * HostArraysPtr.ro_w[i+j*localNx+k*localNx*(def.Ny)] -
                (y1 + abs(y1))/2.0*(-1)* HostArraysPtr.Xi_w[i+(j-1)*localNx+k*localNx*(def.Ny)] * HostArraysPtr.ro_w[i+(j-1)*localNx+k*localNx*(def.Ny)] +
                (y2 - abs(y2))/2.0*(-1)* HostArraysPtr.Xi_w[i+(j+1)*localNx+k*localNx*(def.Ny)] * HostArraysPtr.ro_w[i+(j+1)*localNx+k*localNx*(def.Ny)])/def.hy;


        x2 = -(HostArraysPtr.P_n[i+1+j*localNx+k*localNx*(def.Ny)] - P2)/def.hx;
        x1 = -(P2 - HostArraysPtr.P_n[i-1+j*localNx+k*localNx*(def.Ny)])/def.hx;

        y2 = -(HostArraysPtr.P_n[i+(j+1)*localNx+k*localNx*(def.Ny)] - P2)/def.hy + HostArraysPtr.ro_n[i+j*localNx+k*localNx*(def.Ny)] * (def.g_const);
        y1 = -(P2 - HostArraysPtr.P_n[i+(j-1)*localNx+k*localNx*(def.Ny)])/def.hy + HostArraysPtr.ro_n[i+j*localNx+k*localNx*(def.Ny)] * (def.g_const);

        g1 = (((x2 + abs(x2))/2.0 - (x1 - abs(x1))/2.0)*(-1) * HostArraysPtr.Xi_n[i+j*localNx+k*localNx*(def.Ny)] * HostArraysPtr.ro_n[i+j*localNx+k*localNx*(def.Ny)] -
                (x1 + abs(x1))/2.0*(-1)* HostArraysPtr.Xi_n[i-1+j*localNx+k*localNx*(def.Ny)] * HostArraysPtr.ro_n[i-1+j*localNx+k*localNx*(def.Ny)] +
                (x2 - abs(x2))/2.0*(-1)* HostArraysPtr.Xi_n[i+1+j*localNx+k*localNx*(def.Ny)] * HostArraysPtr.ro_n[i+1+j*localNx+k*localNx*(def.Ny)])/def.hx;

        g2 = (((y2 + abs(y2))/2.0 - (y1 - abs(y1))/2.0)*(-1)* HostArraysPtr.Xi_n[i+j*localNx+k*localNx*(def.Ny)] * HostArraysPtr.ro_n[i+j*localNx+k*localNx*(def.Ny)] -
                (y1 + abs(y1))/2.0*(-1)* HostArraysPtr.Xi_n[i+(j-1)*localNx+k*localNx*(def.Ny)] * HostArraysPtr.ro_n[i+(j-1)*localNx+k*localNx*(def.Ny)] +
                (y2 - abs(y2))/2.0*(-1)* HostArraysPtr.Xi_n[i+(j+1)*localNx+k*localNx*(def.Ny)] * HostArraysPtr.ro_n[i+(j+1)*localNx+k*localNx*(def.Ny)])/def.hy;

		A1 = HostArraysPtr.roS_w[i+j*localNx+k*localNx*(def.Ny)] - (def.dt/def.m[media])*(f1 + f2 + f3);
        A2 = HostArraysPtr.roS_n[i+j*localNx+k*localNx*(def.Ny)] - (def.dt/def.m[media])*(g1 + g2 + g3);

		HostArraysPtr.roS_w_old[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.roS_w[i+j*localNx+k*localNx*(def.Ny)];
		HostArraysPtr.roS_n_old[i+j*localNx+k*localNx*(def.Ny)]= HostArraysPtr.roS_n[i+j*localNx+k*localNx*(def.Ny)];
		HostArraysPtr.roS_w[i+j*localNx+k*localNx*(def.Ny)] = A1;
		HostArraysPtr.roS_n[i+j*localNx+k*localNx*(def.Ny)] = A2;
	}
}

#ifdef THREE_PHASE
void Newton(ptr_Arrays HostArraysPtr, int i, int j, int k, int localNx, consts def)
{
	if ((i!=0) && (i!=localNx-1) && (j!=0) && (j!=(def.Ny)-1) && (((k!=0) && (k!=(def.Nz)-1)) || ((def.Nz)<2)))
	{
		int media = HostArraysPtr.media[i+j*localNx+k*localNx*(def.Ny)];
		double S_w_e, S_g_e, P_k_nw, P_k_gn, AAA, F1, F2, PkS, F1P, F2P, F1S, F2S, det;

		for (int w=1;w<=def.newton_iterations;w++)
		{
			S_w_e = (HostArraysPtr.S_w[i+j*localNx+k*localNx*(def.Ny)] - def.S_wr[media]) / (1. - def.S_wr[media]-def.S_nr[media]-def.S_gr[media]);
			S_g_e = (HostArraysPtr.S_g[i+j*localNx+k*localNx*(def.Ny)] - def.S_gr[media]) / (1. - def.S_wr[media]-def.S_nr[media]-def.S_gr[media]);

			P_k_nw = def.P_d_nw[media]*pow(pow(S_w_e,def.lambda[media]/(1.-def.lambda[media]))-1.,1./def.lambda[media]);
			P_k_gn = def.P_d_nw[media]*pow(pow(1.-S_g_e,def.lambda[media]/(1.-def.lambda[media]))-1.,1./def.lambda[media]);

			AAA = pow(S_e, (((-1.) / def.lambda[media]) - 1.));
			F1 = def.ro0_w * (1. + (def.beta_w) * (HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)] - def.P_atm)) * (1. - HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)]) - HostArraysPtr.roS_w[i+j*localNx+k*localNx*(def.Ny)];
			F2 = def.ro0_n * (1. + (def.beta_n) * (HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)] + P_k - def.P_atm)) * HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)] - HostArraysPtr.roS_n[i+j*localNx+k*localNx*(def.Ny)];
#ifdef TWO_PHASE
			PkS = AAA * (def.P_d[media]) / (def.lambda[media] * (1 - def.S_wr[media]));
#else
			PkS=0;
#endif
			F1P = def.ro0_w * (def.beta_w) * (1 - HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)]);
			F2P = def.ro0_n * (def.beta_n) * HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)];
			F1S = (-1) * (def.ro0_w) * (1 + (def.beta_w) * (HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)] - def.P_atm));
			F2S = def.ro0_n * (1 + (def.beta_n) * (HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)] + P_k - def.P_atm + (HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)] * PkS)));

			det = F1P * F2S - F1S * F2P;

			HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)] - (1 / det) * (F2S * F1 - F1S * F2);
			HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)] - (1 / det) * (F1P * F2 - F2P * F1);
		}  
	}
}
#else
void Newton(ptr_Arrays HostArraysPtr, int i, int j, int k, int localNx, consts def)
{
	if ((i!=0) && (i!=localNx-1) && (j!=0) && (j!=(def.Ny)-1) && (((k!=0) && (k!=(def.Nz)-1)) || ((def.Nz)<2)))
	{
		int media = HostArraysPtr.media[i+j*localNx+k*localNx*(def.Ny)];
		double S_e, P_k, AAA, F1, F2, PkS, F1P, F2P, F1S, F2S, det;

		for (int w=1;w<=def.newton_iterations;w++)
		{
			S_e = (1 - HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)] - def.S_wr[media]) / (1 - def.S_wr[media]);
#ifdef TWO_PHASE
			P_k = def.P_d[media] * pow(S_e, (-1.) / def.lambda[media]);
#else
			P_k=0;
#endif
			AAA = pow(S_e, (((-1.) / def.lambda[media]) - 1.));
			F1 = def.ro0_w * (1. + (def.beta_w) * (HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)] - def.P_atm)) * (1. - HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)]) - HostArraysPtr.roS_w[i+j*localNx+k*localNx*(def.Ny)];
			F2 = def.ro0_n * (1. + (def.beta_n) * (HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)] + P_k - def.P_atm)) * HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)] - HostArraysPtr.roS_n[i+j*localNx+k*localNx*(def.Ny)];
#ifdef TWO_PHASE
			PkS = AAA * (def.P_d[media]) / (def.lambda[media] * (1 - def.S_wr[media]));
#else
			PkS=0;
#endif
			F1P = def.ro0_w * (def.beta_w) * (1 - HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)]);
			F2P = def.ro0_n * (def.beta_n) * HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)];
			F1S = (-1) * (def.ro0_w) * (1 + (def.beta_w) * (HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)] - def.P_atm));
			F2S = def.ro0_n * (1 + (def.beta_n) * (HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)] + P_k - def.P_atm + (HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)] * PkS)));

			det = F1P * F2S - F1S * F2P;

			HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)] - (1 / det) * (F2S * F1 - F1S * F2);
			HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)] - (1 / det) * (F1P * F2 - F2P * F1);
		}  
	}
}
#endif

void Border_Sn(ptr_Arrays HostArraysPtr, int i, int j, int k, int localNx, int rank, int size, consts def)
{
	if ((i == 0) && ((def.Nx)>2))
	{
		HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.S_n[i+1+j*localNx+k*localNx*(def.Ny)];
		return;
	}

	if ((i == localNx - 1) && ((def.Nx)>2))
	{
		HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.S_n[i-1+j*localNx+k*localNx*(def.Ny)];
		return;
	}

	if ((j == (def.Ny) - 1) && ((def.Ny)>2))
	{
		HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.S_n[i+(j-1)*localNx+k*localNx*(def.Ny)];
		return;
	}

	if ((j==0) && ((def.Ny)>2))
	{
		int I=i_to_I(i,rank,size, def);
		if ((I>=(def.Nx)/2-(def.source)) && (I<=(def.Nx)/2+(def.source)) && (k>=(def.Nz)/2-(def.source)) && (k<=(def.Nz)/2+(def.source)))
			HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)] = def.S_n_gr;
		else
			HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)] = 0;
	}

	if ((k == 0) && ((def.Nz)>2))
	{
		HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.S_n[i+j*localNx+(k+1)*localNx*(def.Ny)];
		return;
	}

	if ((k == (def.Nz) - 1) && ((def.Nz)>2))
	{
		HostArraysPtr.S_n[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.S_n[i+j*localNx+(k-1)*localNx*(def.Ny)];
		return;
	}
}

// 
void Border_Pw(ptr_Arrays HostArraysPtr, int i, int j, int k, int localNx, consts def)
{
	if ((i == 0) && ((def.Nx)>2))
	{
		HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.P_w[i+1+j*localNx+k*localNx*(def.Ny)]; 
		return;
	}

	if ((i == localNx - 1) && ((def.Nx)>2))
	{
		HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.P_w[i-1+j*localNx+k*localNx*(def.Ny)];
		return;
	}

	if ((j == (def.Ny) - 1) && ((def.Ny)>2))
	{
		HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.P_w[i+(j-1)*localNx+k*localNx*(def.Ny)] + HostArraysPtr.ro_w[i+localNx*1] * (def.g_const) * (def.hy); ;
		return;
	}

	if ((j==0) && ((def.Ny)>2))
	{
		HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)] = def.P_atm;
		return;
	}

	if ((k == 0) && ((def.Nz)>2))
	{
		HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.P_w[i+j*localNx+(k+1)*localNx*(def.Ny)]; 
		return;
	}

	if ((k == (def.Nz) - 1) && ((def.Nz)>2))
	{
		HostArraysPtr.P_w[i+j*localNx+k*localNx*(def.Ny)] = HostArraysPtr.P_w[i+j*localNx+(k-1)*localNx*(def.Ny)];
		return;
	}
}

// Функция загрузки данных в память хоста
void load_data_to_host(double* HostArrayPtr, double* DevArrayPtr, int localNx, consts def)
{
}

// Функция загрузки данных типа double в память ускорителя
void load_data_to_device(double* HostArrayPtr, double* DevArrayPtr, int localNx, consts def)
{
}

// Функция загрузки данных типа int в память ускорителя
void load_data_to_device_int(int* HostArrayPtr, int* DevArrayPtr, int localNx, consts def)
{
}

// Выделение памяти ускорителя под массив точек расчетной области
void device_memory_allocation(ptr_Arrays* ArraysPtr, double** DevBuffer, int nX, consts def)
{
}

// Освобожение памяти ускорителя из под массива точек расчетной области
void device_memory_free(ptr_Arrays ptDev, double* DevBuffer)
{
}

// Инициализация ускорителя
void device_initialization(int rank, int* blocksX, int* blocksY, int* blocksZ, int localNx, consts def)
{
}

// Загрузка на хост данных для обмена на границе
void load_exchange_data(double* HostArrayPtr, double* DevArrayPtr, double* HostBuffer, double* DevBuffer, int localNx, int blocksY, int blocksZ, int rank, int size, consts def)
{
	for(int j=0;j<(def.Ny);j++)
		for(int k=0;k<(def.Nz);k++)
		{
			HostBuffer[j+(def.Ny)*k]=HostArrayPtr[1+localNx*j+localNx*(def.Ny)*k];
			HostBuffer[j+(def.Ny)*k+(def.Ny)*(def.Nz)]=HostArrayPtr[localNx-2+localNx*j+localNx*(def.Ny)*k];
		}

	/*for(int j=0;j<(def.Ny);j++)
		for(int k=0;k<(def.Nz);k++)
			printf("Buffer j=%d k=%d buffer=%f\n", j, k, HostBuffer[j+(def.Ny)*k]);*/
}

// Загрузка на device данных обмена на границе
void save_exchange_data(double* HostArrayPtr, double* DevArrayPtr, double* HostBuffer, double* DevBuffer, int localNx, int blocksY, int blocksZ, int rank, int size, consts def)
{
	//printf("\nSave\n");
	/*for(int j=0;j<(def.Ny);j++)
		for(int k=0;k<(def.Nz);k++)
			printf("Buffer j=%d k=%d buffer=%f\n", j, k, HostBuffer[j+(def.Ny)*k]);*/

	for(int j=0;j<(def.Ny);j++)
		for(int k=0;k<(def.Nz);k++)
		{
			if(rank!=size-1)
				HostArrayPtr[localNx-1+localNx*j+localNx*(def.Ny)*k]=HostBuffer[j+(def.Ny)*k+(def.Ny)*(def.Nz)];
			if (rank!=0)
				HostArrayPtr[localNx*j+localNx*(def.Ny)*k]=HostBuffer[j+(def.Ny)*k];
		}
}