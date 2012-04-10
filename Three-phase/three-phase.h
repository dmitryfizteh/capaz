#ifndef THREE_PHASE_H
#define THREE_PHASE_H

//  оэффициенты пр€мых, продолжающих функции капилл€рных давлений на границах интервала изменени€ насыщенностей [0,1]
static const double aw[2] = { -396.40, -265.30};
static const double bw[2] = {125.60, 271.10};
static const double ag[2] = {141.70, 353.80};
static const double bg[2] = {1.58, -29.69};
//  лючевые точки интервала изменени€ насыщенностей дл€ вычислени€ проницаемостей и капилл€рных давлений
static const double S_w_range[3] = {0.001, 0.1, 0.99};
static const double S_g_range[3] = {0.001, 0.005, 0.95};

extern void assign_border_P(double* P, double* ro, int i, int j, int k, consts def);

#endif

