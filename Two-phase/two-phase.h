#ifndef TWO_PHASE_H
#define TWO_PHASE_H

extern void Border_Sn(ptr_Arrays HostArraysPtr, int i, int j, int k, int localNx, int rank, int size, consts def);
extern void Border_Pw(ptr_Arrays HostArraysPtr, int i, int j, int k, int localNx, consts def);

extern void data_initialization(ptr_Arrays HostArraysPtr, int* t, int localNx, int localNy, int rank, int size, consts def);
#endif

