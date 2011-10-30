#include "gpu.h"

void gpu_test_nan(double* HostArrayPtr, double* DevArrayPtr, int rank, int size, int localNx, consts def, char* file, int line)
{
#ifdef MY_TEST
	load_data_to_host(HostArrayPtr, DevArrayPtr , localNx, def);
	for(int i=0; i<localNx; i++)
		for(int j=0; j<(def.Ny); j++)
			for(int k=0; k<(def.Nz); k++)
				if(is_active_point(i, localNx, rank, size))
				{
					//std::cout << "i="<<i<<"\tj="<<j<<"\tk="<<k<<"\n";
					test_nan (HostArrayPtr[i+j*localNx+k*localNx*(def.Ny)], file, line);
				}
#endif
}
