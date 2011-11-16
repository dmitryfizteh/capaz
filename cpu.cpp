#include "defines.h"

void ro_P_Xi_calculation(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, consts def, localN locN, int rank, parts_sizes parts, int blocksX, int blocksY, int blocksZ)
{
	for(int i=0;i<(locN.x);i++)
		for(int j=0;j<(locN.y);j++)
			for(int k=0;k<(locN.z);k++)
			{
				if(is_active_point(i, j, k, locN, rank, parts))
				{
					assign_P_Xi(HostArraysPtr,i,j,k,locN,def);
					assign_ro(HostArraysPtr,i,j,k,locN,def);		
				}
			}
}

void u_calculation(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, localN locN, int rank, parts_sizes parts, int blocksX, int blocksY, int blocksZ, consts def)
{
	for(int i=0;i<(locN.x);i++)
		for(int j=0;j<(locN.y);j++)
			for(int k=0;k<(locN.z);k++)
				if(is_active_point(i, j, k, locN, rank, parts))
					assign_u(HostArraysPtr,i,j,k,locN,def);
}

void roS_calculation(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, consts def, double t, localN locN, int rank, parts_sizes parts, int blocksX, int blocksY, int blocksZ)
{
	for(int i=0;i<(locN.x);i++)
		for(int j=0;j<(locN.y);j++)
			for(int k=0;k<(locN.z);k++)
				if(is_active_point(i, j, k, locN, rank, parts))
				{
#ifdef NR
					assign_roS_nr(HostArraysPtr,t,i,j,k,locN,def);
#else
					assign_roS(HostArraysPtr,t,i,j,k,locN,def);
#endif
				}
}

void P_S_calculation(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, consts def, localN locN, int rank, parts_sizes parts, int blocksX, int blocksY, int blocksZ)
{
	for(int i=0;i<(locN.x);i++)
		for(int j=0;j<(locN.y);j++)
			for(int k=0;k<(locN.z);k++)
				if(is_active_point(i, j, k, locN, rank, parts))
					Newton(HostArraysPtr,i,j,k,locN,def);
}

void boundary_conditions(ptr_Arrays HostArraysPtr, ptr_Arrays DevArraysPtr, localN locN, int rank, parts_sizes parts, int blocksX, int blocksY, int blocksZ, consts def)
{
	for(int k = 0; k < (locN.z); k++)
		for(int i = 0; i < (locN.x); i++)
		{
			int j;
			for(j = 1; j < (locN.y); j++)
				if(is_active_point(i, j, k, locN, rank, parts))
				{
					Border_S(HostArraysPtr,i,j,k,locN,rank,parts,def);
					Border_P(HostArraysPtr,i,j,k,locN,def);
				}
			j = 0;
			if(is_active_point(i, j, k, locN, rank, parts))
			{
				Border_S(HostArraysPtr,i,j,k,locN,rank,parts,def);
				Border_P(HostArraysPtr,i,j,k,locN,def);				
			}
		}
}	

void assign_ro(ptr_Arrays HostArraysPtr, int i, int j, int k, localN locN, consts def)
{
	HostArraysPtr.ro_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] = def.ro0_w * (1. + (def.beta_w) * (HostArraysPtr.P_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] - def.P_atm));
	HostArraysPtr.ro_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] = def.ro0_n * (1. + (def.beta_n) * (HostArraysPtr.P_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] - def.P_atm));
#ifdef THREE_PHASE
	HostArraysPtr.ro_g[i+j*(locN.x)+k*(locN.x)*(locN.y)] = def.ro0_g * (1. + (def.beta_g) * (HostArraysPtr.P_g[i+j*(locN.x)+k*(locN.x)*(locN.y)] - def.P_atm));
	test_positive(HostArraysPtr.ro_g[i+j*(locN.x)+k*(locN.x)*(locN.y)], __FILE__, __LINE__);
#endif
	test_positive(HostArraysPtr.ro_w[i+j*(locN.x)+k*(locN.x)*(locN.y)], __FILE__, __LINE__);
	test_positive(HostArraysPtr.ro_n[i+j*(locN.x)+k*(locN.x)*(locN.y)], __FILE__, __LINE__);
}


void assign_u(ptr_Arrays HostArraysPtr, int i, int j, int k, localN locN, consts def)
{
	if ((def.Nx)>2)
	{
		if (i == 0)
		{
#ifdef THREE_PHASE
			HostArraysPtr.ux_g[i+j*(locN.x)+k*(locN.x)*(locN.y)] = HostArraysPtr.Xi_g[i+j*(locN.x)+k*(locN.x)*(locN.y)] * (HostArraysPtr.P_g[i+1+j*(locN.x)+k*(locN.x)*(locN.y)] - HostArraysPtr.P_g[i+j*(locN.x)+k*(locN.x)*(locN.y)]) / (def.hx);
#endif
			HostArraysPtr.ux_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] = HostArraysPtr.Xi_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] * (HostArraysPtr.P_w[i+1+j*(locN.x)+k*(locN.x)*(locN.y)] - HostArraysPtr.P_w[i+j*(locN.x)+k*(locN.x)*(locN.y)]) / (def.hx);
			HostArraysPtr.ux_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] = HostArraysPtr.Xi_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] * (HostArraysPtr.P_n[i+1+j*(locN.x)+k*(locN.x)*(locN.y)] - HostArraysPtr.P_n[i+j*(locN.x)+k*(locN.x)*(locN.y)]) / (def.hx);
		}
		else
		{
			if (i == (locN.x) - 1)
			{
#ifdef THREE_PHASE
				HostArraysPtr.ux_g[i+j*(locN.x)+k*(locN.x)*(locN.y)] = HostArraysPtr.Xi_g[i+j*(locN.x)+k*(locN.x)*(locN.y)] * (HostArraysPtr.P_g[i+j*(locN.x)+k*(locN.x)*(locN.y)] - HostArraysPtr.P_g[i-1+j*(locN.x)+k*(locN.x)*(locN.y)]) / (def.hx);
#endif
				HostArraysPtr.ux_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] = HostArraysPtr.Xi_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] * (HostArraysPtr.P_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] - HostArraysPtr.P_w[i-1+j*(locN.x)+k*(locN.x)*(locN.y)]) / (def.hx);
				HostArraysPtr.ux_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] = HostArraysPtr.Xi_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] * (HostArraysPtr.P_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] - HostArraysPtr.P_n[i-1+j*(locN.x)+k*(locN.x)*(locN.y)]) / (def.hx);
			}
			else
			{
				HostArraysPtr.ux_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] = HostArraysPtr.Xi_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] * ((HostArraysPtr.P_w[i+1+j*(locN.x)+k*(locN.x)*(locN.y)] - HostArraysPtr.P_w[i-1+j*(locN.x)+k*(locN.x)*(locN.y)]) / (2. * (def.hx)));
				HostArraysPtr.ux_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] = HostArraysPtr.Xi_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] * ((HostArraysPtr.P_n[i+1+j*(locN.x)+k*(locN.x)*(locN.y)] - HostArraysPtr.P_n[i-1+j*(locN.x)+k*(locN.x)*(locN.y)]) / (2. * (def.hx)));			
#ifdef THREE_PHASE
				HostArraysPtr.ux_g[i+j*(locN.x)+k*(locN.x)*(locN.y)] = HostArraysPtr.Xi_g[i+j*(locN.x)+k*(locN.x)*(locN.y)] * ((HostArraysPtr.P_g[i+1+j*(locN.x)+k*(locN.x)*(locN.y)] - HostArraysPtr.P_g[i-1+j*(locN.x)+k*(locN.x)*(locN.y)]) / (2. * (def.hx)));
#endif
			}
		}
	}
	else
	{
		HostArraysPtr.ux_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] = 0;
		HostArraysPtr.ux_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] = 0;
#ifdef THREE_PHASE
		HostArraysPtr.ux_g[i+j*(locN.x)+k*(locN.x)*(locN.y)] = 0;
#endif
	}

	if ((def.Ny)>2)
	{
		if (j == 0)
		{
#ifdef THREE_PHASE
			// Если границы области непроницаемые
/*			HostArraysPtr.uy_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] = 0;
			HostArraysPtr.uy_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] = 0;
			HostArraysPtr.uy_g[i+j*(locN.x)+k*(locN.x)*(locN.y)] = 0;
*/			// Если границы области "прозрачные"
			HostArraysPtr.uy_g[i+j*(locN.x)+k*(locN.x)*(locN.y)] = HostArraysPtr.Xi_g[i+j*(locN.x)+k*(locN.x)*(locN.y)] * 
				((HostArraysPtr.P_g[i+(j+1)*(locN.x)+k*(locN.x)*(locN.y)] - HostArraysPtr.ro_g[i+(j+1)*(locN.x)+k*(locN.x)*(locN.y)] * (def.g_const) * (j+1) * (def.hy)) 
				- (HostArraysPtr.P_g[i+j*(locN.x)+k*(locN.x)*(locN.y)] - HostArraysPtr.ro_g[i+j*(locN.x)+k*(locN.x)*(locN.y)] * (def.g_const) * j * (def.hy)))/ def.hy;
			HostArraysPtr.uy_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] = HostArraysPtr.Xi_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] * 
				((HostArraysPtr.P_w[i+(j+1)*(locN.x)+k*(locN.x)*(locN.y)] - HostArraysPtr.ro_w[i+(j+1)*(locN.x)+k*(locN.x)*(locN.y)] * (def.g_const) * (j+1) * (def.hy)) 
				- (HostArraysPtr.P_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] - HostArraysPtr.ro_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] * (def.g_const) * j * (def.hy)))/ def.hy; 
			HostArraysPtr.uy_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] = HostArraysPtr.Xi_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] * 
				((HostArraysPtr.P_n[i+(j+1)*(locN.x)+k*(locN.x)*(locN.y)] - HostArraysPtr.ro_n[i+(j+1)*(locN.x)+k*(locN.x)*(locN.y)] * (def.g_const) * (j+1) * (def.hy)) 
				- (HostArraysPtr.P_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] - HostArraysPtr.ro_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] * (def.g_const) * j * (def.hy)))/ def.hy; 

#else
			HostArraysPtr.uy_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] = HostArraysPtr.Xi_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] * 
				((HostArraysPtr.P_w[i+(j+1)*(locN.x)+k*(locN.x)*(locN.y)] - HostArraysPtr.ro_w[i+(j+1)*(locN.x)+k*(locN.x)*(locN.y)] * (def.g_const) * (j+1) * (def.hy)) 
				- (HostArraysPtr.P_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] - HostArraysPtr.ro_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] * (def.g_const) * j * (def.hy)))/ def.hy; 
			HostArraysPtr.uy_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] = HostArraysPtr.Xi_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] * 
				((HostArraysPtr.P_n[i+(j+1)*(locN.x)+k*(locN.x)*(locN.y)] - HostArraysPtr.ro_n[i+(j+1)*(locN.x)+k*(locN.x)*(locN.y)] * (def.g_const) * (j+1) * (def.hy)) 
				- (HostArraysPtr.P_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] - HostArraysPtr.ro_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] * (def.g_const) * j * (def.hy)))/ def.hy; 
#endif
		}
		else
		{
			if (j == (locN.y) - 1)
			{
#ifdef THREE_PHASE
				HostArraysPtr.uy_g[i+j*(locN.x)+k*(locN.x)*(locN.y)] = HostArraysPtr.Xi_g[i+j*(locN.x)+k*(locN.x)*(locN.y)] * 
					((HostArraysPtr.P_g[i+j*(locN.x)+k*(locN.x)*(locN.y)] - HostArraysPtr.ro_g[i+j*(locN.x)+k*(locN.x)*(locN.y)] * (def.g_const) * j * (def.hy)) 
					- (HostArraysPtr.P_g[i+(j-1)*(locN.x)+k*(locN.x)*(locN.y)] - HostArraysPtr.ro_g[i+(j-1)*(locN.x)+k*(locN.x)*(locN.y)] * (def.g_const) * (j-1) * (def.hy)))/ def.hy;
#endif
				HostArraysPtr.uy_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] = HostArraysPtr.Xi_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] * 
				((HostArraysPtr.P_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] - HostArraysPtr.ro_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] * (def.g_const) * j * (def.hy)) 
					- (HostArraysPtr.P_w[i+(j-1)*(locN.x)+k*(locN.x)*(locN.y)] - HostArraysPtr.ro_w[i+(j-1)*(locN.x)+k*(locN.x)*(locN.y)] * (def.g_const) * (j-1) * (def.hy)))/ def.hy;
				HostArraysPtr.uy_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] = HostArraysPtr.Xi_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] * 
					((HostArraysPtr.P_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] - HostArraysPtr.ro_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] * (def.g_const) * j * (def.hy)) 
					- (HostArraysPtr.P_n[i+(j-1)*(locN.x)+k*(locN.x)*(locN.y)] - HostArraysPtr.ro_n[i+(j-1)*(locN.x)+k*(locN.x)*(locN.y)] * (def.g_const) * (j-1) * (def.hy)))/ def.hy;
			}
			else
			{
				HostArraysPtr.uy_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] = HostArraysPtr.Xi_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] * 
					((HostArraysPtr.P_w[i+(j+1)*(locN.x)+k*(locN.x)*(locN.y)] - HostArraysPtr.ro_w[i+(j+1)*(locN.x)+k*(locN.x)*(locN.y)] * (def.g_const) * (j+1) * (def.hy)) 
					- (HostArraysPtr.P_w[i+(j-1)*(locN.x)+k*(locN.x)*(locN.y)] - HostArraysPtr.ro_w[i+(j-1)*(locN.x)+k*(locN.x)*(locN.y)] * (def.g_const) * (j-1) * (def.hy)))/ (2*def.hy);
				HostArraysPtr.uy_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] = HostArraysPtr.Xi_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] * 
					((HostArraysPtr.P_n[i+(j+1)*(locN.x)+k*(locN.x)*(locN.y)] - HostArraysPtr.ro_n[i+(j+1)*(locN.x)+k*(locN.x)*(locN.y)] * (def.g_const) * (j+1) * (def.hy)) 
					- (HostArraysPtr.P_n[i+(j-1)*(locN.x)+k*(locN.x)*(locN.y)] - HostArraysPtr.ro_n[i+(j-1)*(locN.x)+k*(locN.x)*(locN.y)] * (def.g_const) * (j-1) * (def.hy)))/ (2*def.hy);
#ifdef THREE_PHASE
				HostArraysPtr.uy_g[i+j*(locN.x)+k*(locN.x)*(locN.y)] = HostArraysPtr.Xi_g[i+j*(locN.x)+k*(locN.x)*(locN.y)] * 
					((HostArraysPtr.P_g[i+(j+1)*(locN.x)+k*(locN.x)*(locN.y)] - HostArraysPtr.ro_g[i+(j+1)*(locN.x)+k*(locN.x)*(locN.y)] * (def.g_const) * (j+1) * (def.hy)) 
					- (HostArraysPtr.P_g[i+(j-1)*(locN.x)+k*(locN.x)*(locN.y)] - HostArraysPtr.ro_g[i+(j-1)*(locN.x)+k*(locN.x)*(locN.y)] * (def.g_const) * (j-1) * (def.hy)))/ (2*def.hy);
#endif
			}
		}
	}
	else
	{
		HostArraysPtr.uy_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] = 0;
		HostArraysPtr.uy_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] = 0;
#ifdef THREE_PHASE
		HostArraysPtr.uy_g[i+j*(locN.x)+k*(locN.x)*(locN.y)] = 0;
#endif
	}

	if ((def.Nz)>2)
	{
		if (k == 0)
		{
#ifdef THREE_PHASE
			HostArraysPtr.uz_g[i+j*(locN.x)+k*(locN.x)*(locN.y)] = HostArraysPtr.Xi_g[i+j*(locN.x)+k*(locN.x)*(locN.y)] * (HostArraysPtr.P_g[i+(locN.x)*j+(k+1)*(locN.x)*(locN.y)] - HostArraysPtr.P_g[i+j*(locN.x)+k*(locN.x)*(locN.y)]) / (def.hz);
#endif
			HostArraysPtr.uz_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] = HostArraysPtr.Xi_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] * (HostArraysPtr.P_w[i+(locN.x)*j+(k+1)*(locN.x)*(locN.y)] - HostArraysPtr.P_w[i+j*(locN.x)+k*(locN.x)*(locN.y)]) / (def.hz);
			HostArraysPtr.uz_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] = HostArraysPtr.Xi_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] * (HostArraysPtr.P_n[i+(locN.x)*j+(k+1)*(locN.x)*(locN.y)] - HostArraysPtr.P_n[i+j*(locN.x)+k*(locN.x)*(locN.y)]) / (def.hz);
		}
		else
		{
			if (k == (locN.z) - 1)
			{
#ifdef THREE_PHASE
				HostArraysPtr.uz_g[i+j*(locN.x)+k*(locN.x)*(locN.y)] = HostArraysPtr.Xi_g[i+j*(locN.x)+k*(locN.x)*(locN.y)] * (HostArraysPtr.P_g[i+j*(locN.x)+k*(locN.x)*(locN.y)] - HostArraysPtr.P_g[i+(locN.x)*j+(k-1)*(locN.x)*(locN.y)]) / (def.hz);
#endif
				HostArraysPtr.uz_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] = HostArraysPtr.Xi_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] * (HostArraysPtr.P_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] - HostArraysPtr.P_w[i+(locN.x)*j+(k-1)*(locN.x)*(locN.y)]) / (def.hz);
				HostArraysPtr.uz_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] = HostArraysPtr.Xi_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] * (HostArraysPtr.P_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] - HostArraysPtr.P_n[i+(locN.x)*j+(k-1)*(locN.x)*(locN.y)]) / (def.hz);
			}
			else
			{
				HostArraysPtr.uz_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] = HostArraysPtr.Xi_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] * ((HostArraysPtr.P_w[i+(locN.x)*j+(k+1)*(locN.x)*(locN.y)] - HostArraysPtr.P_w[i+(locN.x)*j+(k-1)*(locN.x)*(locN.y)]) / (2. * (def.hz)));
				HostArraysPtr.uz_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] = HostArraysPtr.Xi_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] * ((HostArraysPtr.P_n[i+(locN.x)*j+(k+1)*(locN.x)*(locN.y)] - HostArraysPtr.P_n[i+(locN.x)*j+(k-1)*(locN.x)*(locN.y)]) / (2. * (def.hz)));
#ifdef THREE_PHASE
				HostArraysPtr.uz_g[i+j*(locN.x)+k*(locN.x)*(locN.y)] = HostArraysPtr.Xi_g[i+j*(locN.x)+k*(locN.x)*(locN.y)] * ((HostArraysPtr.P_g[i+(locN.x)*j+(k+1)*(locN.x)*(locN.y)] - HostArraysPtr.P_g[i+(locN.x)*j+(k-1)*(locN.x)*(locN.y)]) / (2. * (def.hz)));
#endif
			}
		}
	}
	else
	{
		HostArraysPtr.uz_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] = 0;
		HostArraysPtr.uz_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] = 0;
#ifdef THREE_PHASE
		HostArraysPtr.uz_g[i+j*(locN.x)+k*(locN.x)*(locN.y)] = 0;
#endif
	}

	test_nan(HostArraysPtr.ux_w[i+j*(locN.x)+k*(locN.x)*(locN.y)], __FILE__, __LINE__);
	test_nan(HostArraysPtr.ux_n[i+j*(locN.x)+k*(locN.x)*(locN.y)], __FILE__, __LINE__);
	test_nan(HostArraysPtr.uy_w[i+j*(locN.x)+k*(locN.x)*(locN.y)], __FILE__, __LINE__);
	test_nan(HostArraysPtr.uy_n[i+j*(locN.x)+k*(locN.x)*(locN.y)], __FILE__, __LINE__);
	test_nan(HostArraysPtr.uz_w[i+j*(locN.x)+k*(locN.x)*(locN.y)], __FILE__, __LINE__);
	test_nan(HostArraysPtr.uz_n[i+j*(locN.x)+k*(locN.x)*(locN.y)], __FILE__, __LINE__);
#ifdef THREE_PHASE
	test_nan(HostArraysPtr.ux_g[i+j*(locN.x)+k*(locN.x)*(locN.y)], __FILE__, __LINE__);
	test_nan(HostArraysPtr.uy_g[i+j*(locN.x)+k*(locN.x)*(locN.y)], __FILE__, __LINE__);
	test_nan(HostArraysPtr.uz_g[i+j*(locN.x)+k*(locN.x)*(locN.y)], __FILE__, __LINE__);
#endif
}

void assign_roS(ptr_Arrays HostArraysPtr, double t, int i, int j, int k, localN locN, consts def)
{
	if ((i!=0) && (i!=(locN.x)-1) && (j!=0) && (j!=(locN.y)-1) && (((k!=0) && (k!=(locN.z)-1)) || ((def.Nz)<2)))
	{
		int media = HostArraysPtr.media[i+j*(locN.x)+k*(locN.x)*(locN.y)];
		double divgrad1, divgrad2, Tx1, Ty1, Tz1, Tx2, Ty2, Tz2, A1=0, A2=0;

#ifdef THREE_PHASE
		double divgrad3, Tx3, Ty3, Tz3, A3=0;

		double S_n = (1. - HostArraysPtr.S_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] - HostArraysPtr.S_g[i+j*(locN.x)+k*(locN.x)*(locN.y)]);
		double S_w = HostArraysPtr.S_w[i+j*(locN.x)+k*(locN.x)*(locN.y)];
#else
		double S_n = HostArraysPtr.S_n[i+j*(locN.x)+k*(locN.x)*(locN.y)];
		double S_w = 1 - S_n;
#endif
		HostArraysPtr.roS_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] = HostArraysPtr.ro_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] * S_w;
		HostArraysPtr.roS_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] = HostArraysPtr.ro_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] * S_n;
#ifdef THREE_PHASE
		HostArraysPtr.roS_g[i+j*(locN.x)+k*(locN.x)*(locN.y)] = HostArraysPtr.ro_g[i+j*(locN.x)+k*(locN.x)*(locN.y)] * HostArraysPtr.S_g[i+j*(locN.x)+k*(locN.x)*(locN.y)];   
#endif

		if ((def.Nz)<2)
		{
			divgrad1=0;
			divgrad2=0;
			Tz1=0;
			Tz2=0;
#ifdef THREE_PHASE
			divgrad3=0;
			Tz3=0;
#endif
		}
		else
		{
#ifdef THREE_PHASE
			divgrad1 = (def.m[media] * (def.l) * (def.c_w) / 2.) 
				* (HostArraysPtr.ro_w[i+j*(locN.x)+(k+1)*(locN.x)*(locN.y)] * HostArraysPtr.S_w[i+j*(locN.x)+(k+1)*(locN.x)*(locN.y)] 
				- 2 * HostArraysPtr.ro_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] * HostArraysPtr.S_w[i+j*(locN.x)+k*(locN.x)*(locN.y)]
				+ HostArraysPtr.ro_w[i+j*(locN.x)+(k-1)*(locN.x)*(locN.y)] * HostArraysPtr.S_w[i+j*(locN.x)+(k-1)*(locN.x)*(locN.y)]) / ((def.hz) * (def.hz));
			divgrad2 = (def.m[media] * (def.l) * (def.c_n) / 2.) 
				* (HostArraysPtr.ro_n[i+j*(locN.x)+(k+1)*(locN.x)*(locN.y)]
				* (1. - HostArraysPtr.S_w[i+j*(locN.x)+(k+1)*(locN.x)*(locN.y)] - HostArraysPtr.S_g[i+j*(locN.x)+(k+1)*(locN.x)*(locN.y)]) 
				- 2 * HostArraysPtr.ro_n[i+j*(locN.x)+k*(locN.x)*(locN.y)]
				* (1. - HostArraysPtr.S_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] - HostArraysPtr.S_g[i+j*(locN.x)+k*(locN.x)*(locN.y)])  
				+ HostArraysPtr.ro_n[i+j*(locN.x)+(k-1)*(locN.x)*(locN.y)]
				* (1. - HostArraysPtr.S_w[i+j*(locN.x)+(k-1)*(locN.x)*(locN.y)] - HostArraysPtr.S_g[i+j*(locN.x)+(k-1)*(locN.x)*(locN.y)]) ) / ((def.hz) * (def.hz));
			divgrad3 = (def.m[media] * (def.l) * (def.c_g) / 2.) 
				* (HostArraysPtr.ro_g[i+j*(locN.x)+(k+1)*(locN.x)*(locN.y)] * HostArraysPtr.S_g[i+j*(locN.x)+(k+1)*(locN.x)*(locN.y)] 
				- 2 * HostArraysPtr.ro_g[i+j*(locN.x)+k*(locN.x)*(locN.y)] * HostArraysPtr.S_g[i+j*(locN.x)+k*(locN.x)*(locN.y)] 
				+ HostArraysPtr.ro_g[i+j*(locN.x)+(k-1)*(locN.x)*(locN.y)] * HostArraysPtr.S_g[i+j*(locN.x)+(k-1)*(locN.x)*(locN.y)]) / ((def.hz) * (def.hz));
			Tz3 = (HostArraysPtr.ro_g[i+j*(locN.x)+(k+1)*(locN.x)*(locN.y)] * HostArraysPtr.uz_g[i+j*(locN.x)+(k+1)*(locN.x)*(locN.y)] 
				- HostArraysPtr.ro_g[i+j*(locN.x)+(k-1)*(locN.x)*(locN.y)] * HostArraysPtr.uz_g[i+j*(locN.x)+(k-1)*(locN.x)*(locN.y)]) / (2. * (def.hz));            
#else
			divgrad1 = (def.m[media] * (def.l) * (def.c_w) / 2.) 
				* (HostArraysPtr.ro_w[i+j*(locN.x)+(k+1)*(locN.x)*(locN.y)] * (1. - HostArraysPtr.S_n[i+j*(locN.x)+(k+1)*(locN.x)*(locN.y)]) 
				- 2 * HostArraysPtr.ro_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] * (1. - HostArraysPtr.S_n[i+j*(locN.x)+k*(locN.x)*(locN.y)]) 
				+ HostArraysPtr.ro_w[i+j*(locN.x)+(k-1)*(locN.x)*(locN.y)] * (1. - HostArraysPtr.S_n[i+j*(locN.x)+(k-1)*(locN.x)*(locN.y)])) / ((def.hz) * (def.hz));
			divgrad2 = (def.m[media] * (def.l) * (def.c_n) / 2.) 
				* (HostArraysPtr.ro_n[i+j*(locN.x)+(k+1)*(locN.x)*(locN.y)] * HostArraysPtr.S_n[i+j*(locN.x)+(k+1)*(locN.x)*(locN.y)] 
				- 2 * HostArraysPtr.ro_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] * HostArraysPtr.S_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] 
				+ HostArraysPtr.ro_n[i+j*(locN.x)+(k-1)*(locN.x)*(locN.y)] * (HostArraysPtr.S_n[i+j*(locN.x)+(k-1)*(locN.x)*(locN.y)])) / ((def.hz) * (def.hz));
#endif 
			Tz1 = (HostArraysPtr.ro_w[i+j*(locN.x)+(k+1)*(locN.x)*(locN.y)] * HostArraysPtr.uz_w[i+j*(locN.x)+(k+1)*(locN.x)*(locN.y)] 
				- HostArraysPtr.ro_w[i+j*(locN.x)+(k-1)*(locN.x)*(locN.y)] * HostArraysPtr.uz_w[i+j*(locN.x)+(k-1)*(locN.x)*(locN.y)]) / (2. * (def.hz));
			Tz2 = (HostArraysPtr.ro_n[i+j*(locN.x)+(k+1)*(locN.x)*(locN.y)] * HostArraysPtr.uz_n[i+j*(locN.x)+(k+1)*(locN.x)*(locN.y)] 
				- HostArraysPtr.ro_n[i+j*(locN.x)+(k-1)*(locN.x)*(locN.y)] * HostArraysPtr.uz_n[i+j*(locN.x)+(k-1)*(locN.x)*(locN.y)]) / (2. * (def.hz));
		}

#ifdef THREE_PHASE
		divgrad1 += (def.m[media] * (def.l) * (def.c_w) / 2.) *
			((HostArraysPtr.ro_w[i+1+j*(locN.x)+k*(locN.x)*(locN.y)] * (HostArraysPtr.S_w[i+1+j*(locN.x)+k*(locN.x)*(locN.y)])
			- 2 * HostArraysPtr.ro_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] * (HostArraysPtr.S_w[i+j*(locN.x)+k*(locN.x)*(locN.y)]) 
			+ HostArraysPtr.ro_w[i-1+j*(locN.x)+k*(locN.x)*(locN.y)] * (HostArraysPtr.S_w[i-1+j*(locN.x)+k*(locN.x)*(locN.y)])) / ((def.hx) * (def.hx)) +
			(HostArraysPtr.ro_w[i+(j+1)*(locN.x)+k*(locN.x)*(locN.y)] * (HostArraysPtr.S_w[i+(j+1)*(locN.x)+k*(locN.x)*(locN.y)])
			- 2 * HostArraysPtr.ro_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] * (HostArraysPtr.S_w[i+j*(locN.x)+k*(locN.x)*(locN.y)])
			+ HostArraysPtr.ro_w[i+(j-1)*(locN.x)+k*(locN.x)*(locN.y)] * (HostArraysPtr.S_w[i+(j-1)*(locN.x)+k*(locN.x)*(locN.y)])) / ((def.hy) * (def.hy)));

		divgrad2 += (def.m[media] * (def.l) * (def.c_n) / 2.) *
			((HostArraysPtr.ro_n[i+1+j*(locN.x)+k*(locN.x)*(locN.y)] * (1. - HostArraysPtr.S_w[i+1+j*(locN.x)+k*(locN.x)*(locN.y)] - HostArraysPtr.S_g[i+1+j*(locN.x)+k*(locN.x)*(locN.y)])
			- 2 * HostArraysPtr.ro_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] * (1. - HostArraysPtr.S_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] - HostArraysPtr.S_g[i+j*(locN.x)+k*(locN.x)*(locN.y)])
			+ HostArraysPtr.ro_n[i-1+j*(locN.x)+k*(locN.x)*(locN.y)] * (1. - HostArraysPtr.S_w[i-1+j*(locN.x)+k*(locN.x)*(locN.y)] - HostArraysPtr.S_g[i-1+j*(locN.x)+k*(locN.x)*(locN.y)])) / ((def.hx) * (def.hx)) +
			(HostArraysPtr.ro_n[i+(j+1)*(locN.x)+k*(locN.x)*(locN.y)] * (1. - HostArraysPtr.S_w[i+(j+1)*(locN.x)+k*(locN.x)*(locN.y)] - HostArraysPtr.S_g[i+(j+1)*(locN.x)+k*(locN.x)*(locN.y)]) 
			- 2 * HostArraysPtr.ro_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] * (1. - HostArraysPtr.S_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] - HostArraysPtr.S_g[i+j*(locN.x)+k*(locN.x)*(locN.y)]) 
			+ HostArraysPtr.ro_n[i+(j-1)*(locN.x)+k*(locN.x)*(locN.y)] * (1. - HostArraysPtr.S_w[i+(j-1)*(locN.x)+k*(locN.x)*(locN.y)] - HostArraysPtr.S_g[i+(j-1)*(locN.x)+k*(locN.x)*(locN.y)])) / ((def.hy) * (def.hy)));

		divgrad3 += (def.m[media] * (def.l) * (def.c_g) / 2.) *
			((HostArraysPtr.ro_g[i+1+j*(locN.x)+k*(locN.x)*(locN.y)] * (HostArraysPtr.S_g[i+1+j*(locN.x)+k*(locN.x)*(locN.y)]) 
			- 2 * HostArraysPtr.ro_g[i+j*(locN.x)+k*(locN.x)*(locN.y)] * (HostArraysPtr.S_g[i+j*(locN.x)+k*(locN.x)*(locN.y)]) 
			+ HostArraysPtr.ro_g[i-1+j*(locN.x)+k*(locN.x)*(locN.y)] * (HostArraysPtr.S_g[i-1+j*(locN.x)+k*(locN.x)*(locN.y)])) / ((def.hx) * (def.hx)) +
			(HostArraysPtr.ro_g[i+(j+1)*(locN.x)+k*(locN.x)*(locN.y)] * (HostArraysPtr.S_g[i+(j+1)*(locN.x)+k*(locN.x)*(locN.y)]) 
			- 2 * HostArraysPtr.ro_g[i+j*(locN.x)+k*(locN.x)*(locN.y)] * (HostArraysPtr.S_g[i+j*(locN.x)+k*(locN.x)*(locN.y)]) 
			+ HostArraysPtr.ro_g[i+(j-1)*(locN.x)+k*(locN.x)*(locN.y)] * (HostArraysPtr.S_g[i+(j-1)*(locN.x)+k*(locN.x)*(locN.y)])) / ((def.hy) * (def.hy)));

		Tx3 = (HostArraysPtr.ro_g[i+1+j*(locN.x)+k*(locN.x)*(locN.y)] * HostArraysPtr.ux_g[i+1+j*(locN.x)+k*(locN.x)*(locN.y)]
			- HostArraysPtr.ro_g[i-1+j*(locN.x)+k*(locN.x)*(locN.y)] * HostArraysPtr.ux_g[i-1+j*(locN.x)+k*(locN.x)*(locN.y)]) / (2. * (def.hx));
		Ty3 = (HostArraysPtr.ro_g[i+(j+1)*(locN.x)+k*(locN.x)*(locN.y)] * HostArraysPtr.uy_g[i+(j+1)*(locN.x)+k*(locN.x)*(locN.y)] 
			- HostArraysPtr.ro_g[i+(j-1)*(locN.x)+k*(locN.x)*(locN.y)] * HostArraysPtr.uy_g[i+(j-1)*(locN.x)+k*(locN.x)*(locN.y)]) / (2. * (def.hy));
#else
		divgrad1 += (def.m[media] * (def.l) * (def.c_w) / 2.) *
			((HostArraysPtr.ro_w[i+1+j*(locN.x)+k*(locN.x)*(locN.y)] * (1. - HostArraysPtr.S_n[i+1+j*(locN.x)+k*(locN.x)*(locN.y)]) 
			- 2 * HostArraysPtr.ro_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] * (1. - HostArraysPtr.S_n[i+j*(locN.x)+k*(locN.x)*(locN.y)]) 
			+ HostArraysPtr.ro_w[i-1+j*(locN.x)+k*(locN.x)*(locN.y)] * (1. - HostArraysPtr.S_n[i-1+j*(locN.x)+k*(locN.x)*(locN.y)])) / ((def.hx) * (def.hx)) +
			(HostArraysPtr.ro_w[i+(j+1)*(locN.x)+k*(locN.x)*(locN.y)] * (1. - HostArraysPtr.S_n[i+(j+1)*(locN.x)+k*(locN.x)*(locN.y)]) 
			- 2 * HostArraysPtr.ro_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] * (1. - HostArraysPtr.S_n[i+j*(locN.x)+k*(locN.x)*(locN.y)]) 
			+ HostArraysPtr.ro_w[i+(j-1)*(locN.x)+k*(locN.x)*(locN.y)] * (1. - HostArraysPtr.S_n[i+(j-1)*(locN.x)+k*(locN.x)*(locN.y)])) / ((def.hy) * (def.hy)));

		divgrad2 += (def.m[media] * (def.l) * (def.c_n) / 2.) *
			((HostArraysPtr.ro_n[i+1+j*(locN.x)+k*(locN.x)*(locN.y)] * HostArraysPtr.S_n[i+1+j*(locN.x)+k*(locN.x)*(locN.y)] 
			- 2 * HostArraysPtr.ro_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] * HostArraysPtr.S_n[i+j*(locN.x)+k*(locN.x)*(locN.y)]
			+ HostArraysPtr.ro_n[i-1+j*(locN.x)+k*(locN.x)*(locN.y)] * (HostArraysPtr.S_n[i-1+j*(locN.x)+k*(locN.x)*(locN.y)])) / ((def.hx) * (def.hx)) +
			(HostArraysPtr.ro_n[i+(j+1)*(locN.x)+k*(locN.x)*(locN.y)] * HostArraysPtr.S_n[i+(j+1)*(locN.x)+k*(locN.x)*(locN.y)]
			- 2 * HostArraysPtr.ro_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] * HostArraysPtr.S_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] 
			+ HostArraysPtr.ro_n[i+(j-1)*(locN.x)+k*(locN.x)*(locN.y)] * (HostArraysPtr.S_n[i+(j-1)*(locN.x)+k*(locN.x)*(locN.y)])) / ((def.hy) * (def.hy)));
#endif
		Tx1 = (HostArraysPtr.ro_w[i+1+j*(locN.x)+k*(locN.x)*(locN.y)] * HostArraysPtr.ux_w[i+1+j*(locN.x)+k*(locN.x)*(locN.y)] 
			- HostArraysPtr.ro_w[i-1+j*(locN.x)+k*(locN.x)*(locN.y)] * HostArraysPtr.ux_w[i-1+j*(locN.x)+k*(locN.x)*(locN.y)]) / (2. * (def.hx));
		Ty1 = (HostArraysPtr.ro_w[i+(j+1)*(locN.x)+k*(locN.x)*(locN.y)] * HostArraysPtr.uy_w[i+(j+1)*(locN.x)+k*(locN.x)*(locN.y)] 
			- HostArraysPtr.ro_w[i+(j-1)*(locN.x)+k*(locN.x)*(locN.y)] * HostArraysPtr.uy_w[i+(j-1)*(locN.x)+k*(locN.x)*(locN.y)]) / (2. * (def.hy));
		Tx2 = (HostArraysPtr.ro_n[i+1+j*(locN.x)+k*(locN.x)*(locN.y)] * HostArraysPtr.ux_n[i+1+j*(locN.x)+k*(locN.x)*(locN.y)] 
			- HostArraysPtr.ro_n[i-1+j*(locN.x)+k*(locN.x)*(locN.y)] * HostArraysPtr.ux_n[i-1+j*(locN.x)+k*(locN.x)*(locN.y)]) / (2. * (def.hx));
		Ty2 = (HostArraysPtr.ro_n[i+(j+1)*(locN.x)+k*(locN.x)*(locN.y)] * HostArraysPtr.uy_n[i+(j+1)*(locN.x)+k*(locN.x)*(locN.y)]
			- HostArraysPtr.ro_n[i+(j-1)*(locN.x)+k*(locN.x)*(locN.y)] * HostArraysPtr.uy_n[i+(j-1)*(locN.x)+k*(locN.x)*(locN.y)]) / (2. * (def.hy));

			test_arrowhead(Tx1+Ty1+Tz1, divgrad1, __FILE__, __LINE__);
			test_arrowhead(Tx2+Ty2+Tz2, divgrad2, __FILE__, __LINE__);

		if ((t < 2 * (def.dt)) || TWO_LAYERS)
		{
			A1 = HostArraysPtr.roS_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] + ((def.dt) / def.m[media]) * (divgrad1 - Tx1 - Ty1 - Tz1);
			A2 = HostArraysPtr.roS_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] + ((def.dt) / def.m[media]) * (divgrad2 - Tx2 - Ty2 - Tz2);
#ifdef THREE_PHASE
			A3 = HostArraysPtr.roS_g[i+j*(locN.x)+k*(locN.x)*(locN.y)] + ((def.dt) / def.m[media]) * (divgrad3 - Tx3 - Ty3 - Tz3);
#endif
		}
		else
		{
			A1 = (1. / ((def.m[media]) * (def.dt) + 2 * (def.tau))) * (2 * (def.dt) * (def.dt) * (divgrad1 - Tx1 - Ty1 - Tz1) 
				+ ((def.m[media]) * (def.dt) - 2 * (def.tau)) * HostArraysPtr.roS_w_old[i+j*(locN.x)+k*(locN.x)*(locN.y)] 
				+ 4 * (def.tau) * HostArraysPtr.roS_w[i+j*(locN.x)+k*(locN.x)*(locN.y)]);  
			A2 = (1. / ((def.m[media]) * (def.dt) + 2 * (def.tau))) * (2 * (def.dt) * (def.dt) * (divgrad1 - Tx1 - Ty1 - Tz1) 
				+ ((def.m[media]) * (def.dt) - 2 * (def.tau)) * HostArraysPtr.roS_n_old[i+j*(locN.x)+k*(locN.x)*(locN.y)] 
				+ 4 * (def.tau) * HostArraysPtr.roS_n[i+j*(locN.x)+k*(locN.x)*(locN.y)]);  

			test_tau(HostArraysPtr.roS_w_old[i+j*(locN.x)+k*(locN.x)*(locN.y)], HostArraysPtr.roS_w[i+j*(locN.x)+k*(locN.x)*(locN.y)], A1, HostArraysPtr.media[i+j*(locN.x)+k*(locN.x)*(locN.y)], def, __FILE__, __LINE__);
			test_tau(HostArraysPtr.roS_n_old[i+j*(locN.x)+k*(locN.x)*(locN.y)], HostArraysPtr.roS_n[i+j*(locN.x)+k*(locN.x)*(locN.y)], A2, HostArraysPtr.media[i+j*(locN.x)+k*(locN.x)*(locN.y)], def, __FILE__, __LINE__);

#ifdef THREE_PHASE
			A3 = (1. / ((def.m[media]) * (def.dt) + 2 * (def.tau))) * (2 * (def.dt) * (def.dt) * (divgrad1 - Tx1 - Ty1 - Tz1) 
				+ ((def.m[media]) * (def.dt) - 2 * (def.tau)) * HostArraysPtr.roS_g_old[i+j*(locN.x)+k*(locN.x)*(locN.y)] 
				+ 4 * (def.tau) * HostArraysPtr.roS_g[i+j*(locN.x)+k*(locN.x)*(locN.y)]);  
#endif
		}

		HostArraysPtr.roS_w_old[i+j*(locN.x)+k*(locN.x)*(locN.y)] = HostArraysPtr.roS_w[i+j*(locN.x)+k*(locN.x)*(locN.y)];
		HostArraysPtr.roS_n_old[i+j*(locN.x)+k*(locN.x)*(locN.y)] = HostArraysPtr.roS_n[i+j*(locN.x)+k*(locN.x)*(locN.y)];
		HostArraysPtr.roS_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] = A1;
		HostArraysPtr.roS_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] = A2;

		test_positive(HostArraysPtr.roS_w[i+j*(locN.x)+k*(locN.x)*(locN.y)], __FILE__, __LINE__);
		test_positive(HostArraysPtr.roS_n[i+j*(locN.x)+k*(locN.x)*(locN.y)], __FILE__, __LINE__);

#ifdef THREE_PHASE
		HostArraysPtr.roS_g_old[i+j*(locN.x)+k*(locN.x)*(locN.y)] = HostArraysPtr.roS_g[i+j*(locN.x)+k*(locN.x)*(locN.y)];
		HostArraysPtr.roS_g[i+j*(locN.x)+k*(locN.x)*(locN.y)] = A3;
		test_positive(HostArraysPtr.roS_g[i+j*(locN.x)+k*(locN.x)*(locN.y)], __FILE__, __LINE__);
#endif
//		double delta_roS_w = HostArraysPtr.roS_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] - HostArraysPtr.roS_w_old[i+j*(locN.x)+k*(locN.x)*(locN.y)];
//		double delta_roS_n = HostArraysPtr.roS_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] - HostArraysPtr.roS_n_old[i+j*(locN.x)+k*(locN.x)*(locN.y)];
//		double delta_roS_g = HostArraysPtr.roS_g[i+j*(locN.x)+k*(locN.x)*(locN.y)] - HostArraysPtr.roS_g_old[i+j*(locN.x)+k*(locN.x)*(locN.y)];
	}
}

void assign_roS_nr(ptr_Arrays HostArraysPtr, double t, int i, int j, int k, localN locN, consts def)
{
	if ((i!=0) && (i!=(locN.x)-1) && (j!=0) && (j!=(locN.y)-1) && (((k!=0) && (k!=(locN.z)-1)) || ((locN.z)<2)))
	{
		int media = HostArraysPtr.media[i+j*(locN.x)+k*(locN.x)*(locN.y)];

#ifdef THREE_PHASE
		HostArraysPtr.roS_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] = HostArraysPtr.ro_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] 
			* (1. - HostArraysPtr.S_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] - HostArraysPtr.S_g[i+j*(locN.x)+k*(locN.x)*(locN.y)]);
		HostArraysPtr.roS_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] = HostArraysPtr.ro_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] * HostArraysPtr.S_w[i+j*(locN.x)+k*(locN.x)*(locN.y)];
		HostArraysPtr.roS_g[i+j*(locN.x)+k*(locN.x)*(locN.y)] = HostArraysPtr.ro_g[i+j*(locN.x)+k*(locN.x)*(locN.y)] * HostArraysPtr.S_g[i+j*(locN.x)+k*(locN.x)*(locN.y)];

		double Pg = HostArraysPtr.P_g[i+j*(locN.x)+k*(locN.x)*(locN.y)];
		double b1, b2, b3, A3 = 0;
#else
		HostArraysPtr.roS_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] = HostArraysPtr.ro_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] * (1 - HostArraysPtr.S_n[i+j*(locN.x)+k*(locN.x)*(locN.y)]);
		HostArraysPtr.roS_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] = HostArraysPtr.ro_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] * HostArraysPtr.S_n[i+j*(locN.x)+k*(locN.x)*(locN.y)];
#endif
		double Pw = HostArraysPtr.P_w[i+j*(locN.x)+k*(locN.x)*(locN.y)];
		double Pn = HostArraysPtr.P_n[i+j*(locN.x)+k*(locN.x)*(locN.y)];

		double x1, x2, y1, y2, z1, z2, f1, f2, f3, g1, g2, g3, A1=0, A2=0;

		if ((def.Nz) < 2)
		{
			f3 = 0;
			g3 = 0;
#ifdef THREE_PHASE
			b3 = 0;
#endif
		}
		else
		{
			z2 = -(HostArraysPtr.P_w[i+j*(locN.x)+(k+1)*(locN.x)*(locN.y)] - Pw)/def.hz;
			z1 = -(Pw - HostArraysPtr.P_w[i+j*(locN.x)+(k-1)*(locN.x)*(locN.y)])/def.hz;

			f3 = (((z2 + abs(z2))/2.0 - (z1 - abs(z1))/2.0)*(-1) * HostArraysPtr.Xi_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] * HostArraysPtr.ro_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] -
                      (z1 + abs(z1))/2.0*(-1)* HostArraysPtr.Xi_w[i+j*(locN.x)+(k-1)*(locN.x)*(locN.y)] * HostArraysPtr.ro_w[i+j*(locN.x)+(k-1)*(locN.x)*(locN.y)] +
                      (z2 - abs(z2))/2.0*(-1)* HostArraysPtr.Xi_w[i+j*(locN.x)+(k+1)*(locN.x)*(locN.y)] * HostArraysPtr.ro_w[i+j*(locN.x)+(k+1)*(locN.x)*(locN.y)])/def.hz;

			z2 = -(HostArraysPtr.P_n[i+j*(locN.x)+(k+1)*(locN.x)*(locN.y)] - Pn)/def.hz;
			z1 = -(Pn - HostArraysPtr.P_n[i+j*(locN.x)+(k-1)*(locN.x)*(locN.y)])/def.hz;

			g3 = (((z2 + abs(z2))/2.0 - (z1 - abs(z1))/2.0)*(-1) * HostArraysPtr.Xi_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] * HostArraysPtr.ro_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] -
                      (z1 + abs(z1))/2.0*(-1)* HostArraysPtr.Xi_n[i+j*(locN.x)+(k-1)*(locN.x)*(locN.y)] * HostArraysPtr.ro_n[i+j*(locN.x)+(k-1)*(locN.x)*(locN.y)] +
                      (z2 - abs(z2))/2.0*(-1)* HostArraysPtr.Xi_n[i+j*(locN.x)+(k+1)*(locN.x)*(locN.y)] * HostArraysPtr.ro_n[i+j*(locN.x)+(k+1)*(locN.x)*(locN.y)])/def.hz;
#ifdef THREE_PHASE
			z2 = -(HostArraysPtr.P_g[i+j*(locN.x)+(k+1)*(locN.x)*(locN.y)] - Pg)/def.hz;
			z1 = -(Pg - HostArraysPtr.P_g[i+j*(locN.x)+(k-1)*(locN.x)*(locN.y)])/def.hz;

			b3 = (((z2 + abs(z2))/2.0 - (z1 - abs(z1))/2.0)*(-1) * HostArraysPtr.Xi_g[i+j*(locN.x)+k*(locN.x)*(locN.y)] * HostArraysPtr.ro_g[i+j*(locN.x)+k*(locN.x)*(locN.y)] -
				(z1 + abs(z1))/2.0*(-1)* HostArraysPtr.Xi_g[i+j*(locN.x)+(k-1)*(locN.x)*(locN.y)] * HostArraysPtr.ro_g[i+j*(locN.x)+(k-1)*(locN.x)*(locN.y)] +
				(z2 - abs(z2))/2.0*(-1)* HostArraysPtr.Xi_g[i+j*(locN.x)+(k+1)*(locN.x)*(locN.y)] * HostArraysPtr.ro_g[i+j*(locN.x)+(k+1)*(locN.x)*(locN.y)])/def.hz;
#endif
		}

		x2 = -(HostArraysPtr.P_w[i+1+j*(locN.x)+k*(locN.x)*(locN.y)] - Pw)/def.hx;
        x1 = -(Pw - HostArraysPtr.P_w[i-1+j*(locN.x)+k*(locN.x)*(locN.y)])/def.hx;

        y2 = -(HostArraysPtr.P_w[i+(j+1)*(locN.x)+k*(locN.x)*(locN.y)] - Pw)/def.hy + HostArraysPtr.ro_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] * (def.g_const);
        y1 = -(Pw - HostArraysPtr.P_w[i+(j-1)*(locN.x)+k*(locN.x)*(locN.y)])/def.hy + HostArraysPtr.ro_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] * (def.g_const);

        f1 = (((x2 + abs(x2))/2.0 - (x1 - abs(x1))/2.0)*(-1) * HostArraysPtr.Xi_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] * HostArraysPtr.ro_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] -
                (x1 + abs(x1))/2.0*(-1)* HostArraysPtr.Xi_w[i-1+j*(locN.x)+k*(locN.x)*(locN.y)] * HostArraysPtr.ro_w[i-1+j*(locN.x)+k*(locN.x)*(locN.y)] +
                (x2 - abs(x2))/2.0*(-1)* HostArraysPtr.Xi_w[i+1+j*(locN.x)+k*(locN.x)*(locN.y)] * HostArraysPtr.ro_w[i+1+j*(locN.x)+k*(locN.x)*(locN.y)])/def.hx;

        f2 = (((y2 + abs(y2))/2.0 - (y1 - abs(y1))/2.0)*(-1)* HostArraysPtr.Xi_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] * HostArraysPtr.ro_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] -
                (y1 + abs(y1))/2.0*(-1)* HostArraysPtr.Xi_w[i+(j-1)*(locN.x)+k*(locN.x)*(locN.y)] * HostArraysPtr.ro_w[i+(j-1)*(locN.x)+k*(locN.x)*(locN.y)] +
                (y2 - abs(y2))/2.0*(-1)* HostArraysPtr.Xi_w[i+(j+1)*(locN.x)+k*(locN.x)*(locN.y)] * HostArraysPtr.ro_w[i+(j+1)*(locN.x)+k*(locN.x)*(locN.y)])/def.hy;


        x2 = -(HostArraysPtr.P_n[i+1+j*(locN.x)+k*(locN.x)*(locN.y)] - Pn)/def.hx;
        x1 = -(Pn - HostArraysPtr.P_n[i-1+j*(locN.x)+k*(locN.x)*(locN.y)])/def.hx;

        y2 = -(HostArraysPtr.P_n[i+(j+1)*(locN.x)+k*(locN.x)*(locN.y)] - Pn)/def.hy + HostArraysPtr.ro_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] * (def.g_const);
        y1 = -(Pn - HostArraysPtr.P_n[i+(j-1)*(locN.x)+k*(locN.x)*(locN.y)])/def.hy + HostArraysPtr.ro_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] * (def.g_const);

        g1 = (((x2 + abs(x2))/2.0 - (x1 - abs(x1))/2.0)*(-1) * HostArraysPtr.Xi_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] * HostArraysPtr.ro_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] -
                (x1 + abs(x1))/2.0*(-1)* HostArraysPtr.Xi_n[i-1+j*(locN.x)+k*(locN.x)*(locN.y)] * HostArraysPtr.ro_n[i-1+j*(locN.x)+k*(locN.x)*(locN.y)] +
                (x2 - abs(x2))/2.0*(-1)* HostArraysPtr.Xi_n[i+1+j*(locN.x)+k*(locN.x)*(locN.y)] * HostArraysPtr.ro_n[i+1+j*(locN.x)+k*(locN.x)*(locN.y)])/def.hx;

        g2 = (((y2 + abs(y2))/2.0 - (y1 - abs(y1))/2.0)*(-1)* HostArraysPtr.Xi_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] * HostArraysPtr.ro_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] -
                (y1 + abs(y1))/2.0*(-1)* HostArraysPtr.Xi_n[i+(j-1)*(locN.x)+k*(locN.x)*(locN.y)] * HostArraysPtr.ro_n[i+(j-1)*(locN.x)+k*(locN.x)*(locN.y)] +
                (y2 - abs(y2))/2.0*(-1)* HostArraysPtr.Xi_n[i+(j+1)*(locN.x)+k*(locN.x)*(locN.y)] * HostArraysPtr.ro_n[i+(j+1)*(locN.x)+k*(locN.x)*(locN.y)])/def.hy;

		A1 = HostArraysPtr.roS_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] - (def.dt/def.m[media])*(f1 + f2 + f3);
        A2 = HostArraysPtr.roS_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] - (def.dt/def.m[media])*(g1 + g2 + g3);

		HostArraysPtr.roS_w_old[i+j*(locN.x)+k*(locN.x)*(locN.y)] = HostArraysPtr.roS_w[i+j*(locN.x)+k*(locN.x)*(locN.y)];
		HostArraysPtr.roS_n_old[i+j*(locN.x)+k*(locN.x)*(locN.y)]= HostArraysPtr.roS_n[i+j*(locN.x)+k*(locN.x)*(locN.y)];
		HostArraysPtr.roS_w[i+j*(locN.x)+k*(locN.x)*(locN.y)] = A1;
		HostArraysPtr.roS_n[i+j*(locN.x)+k*(locN.x)*(locN.y)] = A2;

		test_positive(HostArraysPtr.roS_w[i+j*(locN.x)+k*(locN.x)*(locN.y)], __FILE__, __LINE__);
		test_positive(HostArraysPtr.roS_n[i+j*(locN.x)+k*(locN.x)*(locN.y)], __FILE__, __LINE__);

#ifdef THREE_PHASE
		x2 = -(HostArraysPtr.P_g[i+1+j*(locN.x)+k*(locN.x)*(locN.y)] - Pg)/def.hx;
		x1 = -(Pg - HostArraysPtr.P_g[i-1+j*(locN.x)+k*(locN.x)*(locN.y)])/def.hx;

		y2 = -(HostArraysPtr.P_g[i+(j+1)*(locN.x)+k*(locN.x)*(locN.y)] - Pg)/def.hy + HostArraysPtr.ro_g[i+j*(locN.x)+k*(locN.x)*(locN.y)] * (def.g_const);
		y1 = -(Pg - HostArraysPtr.P_g[i+(j-1)*(locN.x)+k*(locN.x)*(locN.y)])/def.hy + HostArraysPtr.ro_g[i+j*(locN.x)+k*(locN.x)*(locN.y)] * (def.g_const);

		b1 = (((x2 + abs(x2))/2.0 - (x1 - abs(x1))/2.0)*(-1) * HostArraysPtr.Xi_g[i+j*(locN.x)+k*(locN.x)*(locN.y)] * HostArraysPtr.ro_g[i+j*(locN.x)+k*(locN.x)*(locN.y)] -
			(x1 + abs(x1))/2.0*(-1)* HostArraysPtr.Xi_g[i-1+j*(locN.x)+k*(locN.x)*(locN.y)] * HostArraysPtr.ro_g[i-1+j*(locN.x)+k*(locN.x)*(locN.y)] +
			(x2 - abs(x2))/2.0*(-1)* HostArraysPtr.Xi_g[i+1+j*(locN.x)+k*(locN.x)*(locN.y)] * HostArraysPtr.ro_g[i+1+j*(locN.x)+k*(locN.x)*(locN.y)])/def.hx;

		b2 = (((y2 + abs(y2))/2.0 - (y1 - abs(y1))/2.0)*(-1)* HostArraysPtr.Xi_g[i+j*(locN.x)+k*(locN.x)*(locN.y)] * HostArraysPtr.ro_g[i+j*(locN.x)+k*(locN.x)*(locN.y)] -
			(y1 + abs(y1))/2.0*(-1)* HostArraysPtr.Xi_g[i+(j-1)*(locN.x)+k*(locN.x)*(locN.y)] * HostArraysPtr.ro_g[i+(j-1)*(locN.x)+k*(locN.x)*(locN.y)] +
			(y2 - abs(y2))/2.0*(-1)* HostArraysPtr.Xi_g[i+(j+1)*(locN.x)+k*(locN.x)*(locN.y)] * HostArraysPtr.ro_g[i+(j+1)*(locN.x)+k*(locN.x)*(locN.y)])/def.hy;

		A3 = HostArraysPtr.roS_g[i+j*(locN.x)+k*(locN.x)*(locN.y)] - (def.dt/def.m[media])*(b1 + b2 + b3);

		HostArraysPtr.roS_g_old[i+j*(locN.x)+k*(locN.x)*(locN.y)]= HostArraysPtr.roS_g[i+j*(locN.x)+k*(locN.x)*(locN.y)];
		HostArraysPtr.roS_g[i+j*(locN.x)+k*(locN.x)*(locN.y)] = A3;

		test_positive(HostArraysPtr.roS_g[i+j*(locN.x)+k*(locN.x)*(locN.y)], __FILE__, __LINE__);
#endif
	}
}

// Функция загрузки данных в память хоста
void load_data_to_host(double* HostArrayPtr, double* DevArrayPtr, localN locN, consts def)
{
}

// Функция загрузки данных типа double в память ускорителя
void load_data_to_device(double* HostArrayPtr, double* DevArrayPtr, localN locN, consts def)
{
}

// Функция загрузки данных типа int в память ускорителя
void load_data_to_device_int(int* HostArrayPtr, int* DevArrayPtr, localN locN, consts def)
{
}

// Выделение памяти ускорителя под массив точек расчетной области
void device_memory_allocation(ptr_Arrays* ArraysPtr, double** DevBuffer, localN locN, consts def)
{
}

// Освобожение памяти ускорителя из под массива точек расчетной области
void device_memory_free(ptr_Arrays ptDev, double* DevBuffer)
{
}

// Инициализация ускорителя
void device_initialization(int rank, int* blocksX, int* blocksY, int* blocksZ, localN locN, consts def)
{
}

// Финализация ускорителя
void device__finalization(void)
{
}

// Загрузка на хост данных для обмена на границе
void load_exchange_data(double* HostArrayPtr, double* DevArrayPtr, double* HostBuffer, double* DevBuffer, localN locN, int blocksY, int blocksZ, int rank, parts_sizes parts, consts def)
{
	for(int j=0;j<(locN.y);j++)
		for(int k=0;k<(locN.z);k++)
		{
			HostBuffer[j+(locN.y)*k]=HostArrayPtr[1+(locN.x)*j+(locN.x)*(locN.y)*k];
			HostBuffer[j+(locN.y)*k+(locN.y)*(locN.z)]=HostArrayPtr[(locN.x)-2+(locN.x)*j+(locN.x)*(locN.y)*k];

			test_nan(HostBuffer[j+(locN.y)*k], __FILE__, __LINE__);
			test_nan(HostBuffer[j+(locN.y)*k+(locN.y)*(locN.z)], __FILE__, __LINE__);
		}

	/*for(int j=0;j<(locN.y);j++)
		for(int k=0;k<(locN.z);k++)
			printf("Buffer j=%d k=%d buffer=%f\n", j, k, HostBuffer[j+(locN.y)*k]);*/
}

// Загрузка на device данных обмена на границе
void save_exchange_data(double* HostArrayPtr, double* DevArrayPtr, double* HostBuffer, double* DevBuffer, localN locN, int blocksY, int blocksZ, int rank, parts_sizes parts, consts def)
{
	//printf("\nSave\n");
	/*for(int j=0;j<(locN.y);j++)
		for(int k=0;k<(locN.z);k++)
			printf("Buffer j=%d k=%d buffer=%f\n", j, k, HostBuffer[j+(locN.y)*k]);*/

	for(int j=0;j<(locN.y);j++)
		for(int k=0;k<(locN.z);k++)
		{
			if(rank!=parts.x-1)
			{
				HostArrayPtr[(locN.x)-1+(locN.x)*j+(locN.x)*(locN.y)*k]=HostBuffer[j+(locN.y)*k+(locN.y)*(locN.z)];
				test_nan(HostArrayPtr[(locN.x)-1+(locN.x)*j+(locN.x)*(locN.y)*k], __FILE__, __LINE__);
			}
			if (rank!=0)
			{
				HostArrayPtr[(locN.x)*j+(locN.x)*(locN.y)*k]=HostBuffer[j+(locN.y)*k];
				test_nan(HostArrayPtr[(locN.x)*j+(locN.x)*(locN.y)*k], __FILE__, __LINE__);
			}
		}
}
