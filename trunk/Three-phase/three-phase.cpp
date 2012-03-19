#include "../defines.h"
#include "three-phase.h"

// Коэффициенты прямых, продолжающих функции капиллярных давлений на границах интервала изменения насыщенностей [0,1]
static const double aw[2] = { -396.40, -265.30};
static const double bw[2] = {125.60, 271.10};
static const double ag[2] = {141.70, 353.80};
static const double bg[2] = {1.58, -29.69};
// Ключевые точки интервала изменения насыщенностей для вычисления проницаемостей и капиллярных давлений
static const double S_w_range[3] = {0.001, 0.1, 0.99};
static const double S_g_range[3] = {0.001, 0.005, 0.95};

// Функции вычисления капиллярных давлений
// По краям интервала [0, 1] функции капиллярных давлений гладко заменяем линейными, производные меняются соответственно.
// Описание можно посмотреть в файле mathcad.
static double assign_P_k_nw(double S_w_e, int media, consts def)
{
	double P_k_nw = 0;
	double A = def.lambda[media];

	if (S_w_e <= S_w_range[1])
	{
		P_k_nw = aw[0] * S_w_e + bw[0];
	}
	else if (S_w_e >= S_w_range[2])
	{
		P_k_nw = aw[1] * S_w_e + bw[1];
	}
	else
	{
		P_k_nw = def.P_d_nw[media] * pow((pow(S_w_e, A / (1. - A)) - 1.), 1. / A);
	}

	return P_k_nw;
}

static double assign_P_k_gn(double S_g_e, int media, consts def)
{
	double P_k_gn = 0;
	double A = def.lambda[media];

	if (S_g_e <= S_g_range[1])
	{
		P_k_gn = ag[0] * S_g_e + bg[0];
	}
	else if (S_g_e >= S_g_range[2])
	{
		P_k_gn = ag[1] * S_g_e + bg[1];
	}
	else
	{
		P_k_gn = def.P_d_gn[media] * pow(pow((1. - S_g_e), A / (1. - A)) - 1., 1. / A);
	}

	return P_k_gn;
}

// Функции вычисления производных капиллярных давлений по насыщенностям
static double assign_P_k_nw_S(double S_w_e, int media, consts def)
{
	double PkSw = 0;
	double A = def.lambda[media];                                                                                                                                                                                                                                                                  /*2*/

	if (S_w_e <= S_w_range[1])
	{
		PkSw = aw[0];
	}
	else if (S_w_e >= S_w_range[2])
	{
		PkSw = aw[1];
	}
	else
	{
		PkSw = def.P_d_nw[media] * pow(pow(S_w_e, A / (1. - A)) - 1., 1. / A - 1.) * pow(S_w_e, (A / (1. - A) - 1.)) / (1. - A)
			/ (1. - def.S_wr[media] - def.S_nr[media] - def.S_gr[media]);
	}

	return PkSw;
}

static double assign_P_k_gn_S(double S_g_e, int media, consts def)
{
	double PkSn = 0;
	double A = def.lambda[media];

	if (S_g_e <= S_g_range[1])
	{
		PkSn = (-1) * ag[0];
	}
	else if (S_g_e >= S_g_range[2])
	{
		PkSn = (-1) * ag[1];
	}
	else
	{
		PkSn = def.P_d_gn[media] * pow(pow(1. - S_g_e, A / (1. - A)) - 1., 1. / A - 1.) * pow(1. - S_g_e, A / (1. - A) - 1.) / (1. - A)
			/ (1. - def.S_wr[media] - def.S_nr[media] - def.S_gr[media]);
	}
	return PkSn;
}

// Функции вычисления относительных проницаемостей
static double assign_k_w(double S_w_e, int media, consts def)
{
	double A = def.lambda[media];
	double k_w = 0;

	if (S_w_e >= S_w_range[0])
	{
		k_w = pow(S_w_e, 0.5) * pow(1. - pow(1. - pow(S_w_e, A / (A - 1.)), (A - 1.) / A), 2.);
	}

	return k_w;
}

static double assign_k_g(double S_g_e, int media, consts def)
{
	double A = def.lambda[media];
	double k_g = 0;

	if (S_g_e >= S_g_range[0])
	{
		k_g = pow(S_g_e, 0.5) * pow(1. - pow(1. - S_g_e, A / (A - 1.)), 2. * (A - 1.) / A);
	}

	return k_g;
}

static double assign_k_n(double S_w_e, double S_n_e, int media, consts def)
{
	double A = def.lambda[media];
	double k_n = 0;
	double S_g_e = 1. - S_w_e - S_n_e;

	if (S_n_e >= 1e-3)
	{
		double k_n_w = pow(1. - S_w_e, 0.5) * pow(1. - pow(S_w_e, A / (A - 1.)), 2. * (A - 1.) / A);
		double k_n_g = pow(S_n_e, 0.5) * pow(1. - pow(1. - pow(S_n_e, A / (A - 1.)), (A - 1.) / A), 2.);
		k_n = S_n_e * k_n_w * k_n_g / (1 - S_w_e) / (1 - S_g_e);
	}

	return k_n;
}

//Функция вычисления значений давлений, плотностей и коэффициентов в законе Дарси в точке (i,j,k) среды media,
//исходя из известных значений основных параметров (Pw,Sw,Sn)
//1. Запоминаем, с какой именно из сред работаем
//2. Вычисление значения насыщенности фазы n из условия равенства насыщенностей в сумме единице
//3. Вычисление эффективных насыщенностей по формулам модели трехфазной фильтрации
//4. Вычисление относительных фазовых проницаемостей в соответствии с приближением Стоуна в модификации Азиза и Сеттари
//5. Вычисление капиллярных давлений в соответствии с приближенной моделью Паркера
//6. Вычисление фазовых давлений c помощью капиллярных
//7. Вычисление коэффициентов закона Дарси

void assign_P_Xi(ptr_Arrays HostArraysPtr, int i, int j, int k, consts def)
{

	int media = HostArraysPtr.media[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)];
	double k_w, k_g, k_n, P_k_nw, P_k_gn;
	double A = def.lambda[media];
	double S_w_e = (HostArraysPtr.S_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] - def.S_wr[media]) / (1. - def.S_wr[media] - def.S_nr[media] - def.S_gr[media]);
	double S_n_e = (HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] - def.S_nr[media]) / (1. - def.S_nr[media] - def.S_nr[media] - def.S_nr[media]);
	double S_g_e = 1. - S_w_e - S_n_e;

	k_w = assign_k_w(S_w_e, media, def);
	k_g = assign_k_g(S_g_e, media, def);
	k_n = assign_k_n(S_w_e, S_n_e, media, def);

	P_k_nw = assign_P_k_nw(S_w_e, media, def);
	P_k_gn = assign_P_k_gn(S_g_e, media, def);

	HostArraysPtr.P_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] + P_k_nw;
	HostArraysPtr.P_g[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = HostArraysPtr.P_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] + P_k_gn;

	HostArraysPtr.Xi_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = (-1.) * (def.K[media]) * k_w / def.mu_w;
	HostArraysPtr.Xi_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = (-1.) * (def.K[media]) * k_n / def.mu_n;
	HostArraysPtr.Xi_g[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = (-1.) * (def.K[media]) * k_g / def.mu_g;

	test_positive(HostArraysPtr.P_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)], __FILE__, __LINE__);
	test_positive(HostArraysPtr.P_g[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)], __FILE__, __LINE__);
	test_nan(HostArraysPtr.Xi_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)], __FILE__, __LINE__);
	test_nan(HostArraysPtr.Xi_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)], __FILE__, __LINE__);
	test_nan(HostArraysPtr.Xi_g[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)], __FILE__, __LINE__);
}

// Вспомогательная функции для метода Ньютона:
// Нахождение обратной матрицы 3*3;
static void reverse_matrix(double* a)
{
	int n = 3;
	double b[9], det = 0;

	// Вычисление дополнительных миноров матрицы
	for(int j = 0; j < n; j++)
		for(int i = 0; i < n; i++)
		{
			b[i + n * j] = a[(i + 1) % n + n * ((j + 1) % n)] * a[(i + 2) % n + n * ((j + 2) % n)]
			- a[(i + 2) % n + n * ((j + 1) % n)] * a[(i + 1) % n + n * ((j + 2) % n)];
		}

		// Нахождение детерминанта матрицы 3*3;
		for(int i = 0; i < n; i++)
		{
			det += a[i] * b[i];
		}
		test_nan(det, __FILE__, __LINE__);

		// Транспонирование и деление на детерминант
		for(int j = 0; j < n; j++)
			for(int i = 0; i < n; i++)
			{
				a[i + n * j] = b[j + n * i] / det;
				test_nan(a[i + n * j], __FILE__, __LINE__);
			}
}

//Функция решения системы 3*3 на основные параметры (Pn,Sw,Sg) методом Ньютона в точке (i,j,k) среды media
//1. Вычисление эффективных насыщенностей
//2. Переобозначение степенного коэффициента
//3. Вычисление капиллярных давлений
//4. Вычисление насыщенности фазы n
//5. Нахождение значения трех функций системы
//6. Вычисление частных производных производных капиллярных давлений по насыщенностям
//7. Вычисление матрицы частных производных
//8. Вычисление детерминанта матрицы частных производных
//9. Получение решения системы методом Крамера в явном виде
void Newton(ptr_Arrays HostArraysPtr, int i, int j, int k, consts def)
{
	if ((i != 0) && (i != (def.locNx) - 1) && (j != 0) && (j != (def.locNy) - 1) && (((k != 0) && (k != (def.locNz) - 1)) || ((def.locNz) < 2)))
	{
		int media = HostArraysPtr.media[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)];
		double S_w_e, S_g_e, S_n_e, P_k_nw, P_k_gn, PkSw, PkSn, Sg, F1, F2, F3;
		double dF[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};

		for (int w = 1; w <= def.newton_iterations; w++)
		{
			S_w_e = (HostArraysPtr.S_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] - def.S_wr[media]) / (1. - def.S_wr[media] - def.S_nr[media] - def.S_gr[media]);
			S_n_e = (HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] - def.S_nr[media]) / (1. - def.S_nr[media] - def.S_nr[media] - def.S_nr[media]);
			S_g_e = 1. - S_w_e - S_n_e;

			P_k_nw = assign_P_k_nw(S_w_e, media, def);
			P_k_gn = assign_P_k_gn(S_g_e, media, def);
			PkSw = assign_P_k_nw_S(S_w_e, media, def);
			PkSn = assign_P_k_gn_S(S_g_e, media, def);

			Sg = 1. - HostArraysPtr.S_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] - HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)];

			F1 = def.ro0_w * (1. + (def.beta_w) * (HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] - def.P_atm))
			     * HostArraysPtr.S_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] - HostArraysPtr.roS_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)];
			F2 = def.ro0_n * (1. + (def.beta_n) * (HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] + P_k_nw - def.P_atm))
			     * HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] - HostArraysPtr.roS_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)];
			F3 = def.ro0_g * (HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] + P_k_nw + P_k_gn) / def.P_atm
			     * Sg - HostArraysPtr.roS_g[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)];

			// Матрица частных производных. Индексу от 0 до 8 соответствуют F1P, F1Sw, F1Sn, F2P, F2Sw, F2Sn, F3P, F3Sw, F3Sn
			dF[0] = def.ro0_w * def.beta_w * HostArraysPtr.S_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)];
			dF[3] = def.ro0_n * def.beta_n * HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)];
			dF[6] = def.ro0_g * Sg / def.P_atm;
			dF[1] = def.ro0_w * (1 + def.beta_w * (HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] - def.P_atm));
			dF[4] = def.ro0_n * (1. + (def.beta_n) * PkSw) * HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)];
			dF[7] = (-1) * def.ro0_g * (HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] + P_k_nw + P_k_gn - Sg * (PkSn + PkSw)) / def.P_atm;
			dF[2] = 0;
			dF[5] = def.ro0_n * (1. + def.beta_n * (HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] + P_k_nw - def.P_atm));
			dF[8] = (-1) * def.ro0_g * (HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] + P_k_nw + P_k_gn - Sg * PkSn) / def.P_atm;

			reverse_matrix(dF);

			HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)]
			        - (dF[0] * F1 + dF[1] * F2 + dF[2] * F3);
			HostArraysPtr.S_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = HostArraysPtr.S_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)]
			        - (dF[3] * F1 + dF[4] * F2 + dF[5] * F3);
			HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)]
			        - (dF[6] * F1 + dF[7] * F2 + dF[8] * F3);
		}

		test_S(HostArraysPtr.S_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)], __FILE__, __LINE__);
		test_S(HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)], __FILE__, __LINE__);
		test_positive(HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)], __FILE__, __LINE__);
	}
}

//Задание граничных условий отдельно для (Sw,Sg),Pn

// Задание граничных условий с меньшим числом проверок, но с введением дополнительных переменных
void Border_S(ptr_Arrays HostArraysPtr, int i, int j, int k, consts def)
{
	int i1 = i, j1 = j, k1 = k;

	if (i == 0)
	{
		i1 ++;
	}
	if (i == (def.locNx) - 1)
	{
		i1 --;
	}
	if (j == 0)
	{
		j1 ++;
	}
	if (j == (def.locNy) - 1)
	{
		j1 --;
	}
	if ((k == 0) && ((def.locNz) > 2))
	{
		k1 ++;
	}
	if ((k == (def.locNz) - 1) && ((def.locNz) > 2))
	{
		k1 --;
	}

	if ((j != 0) || ((def.source) <= 0))
	{
		HostArraysPtr.S_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = HostArraysPtr.S_w[i1 + j1 * (def.locNx) + k1 * (def.locNx) * (def.locNy)];
		HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = HostArraysPtr.S_n[i1 + j1 * (def.locNx) + k1 * (def.locNx) * (def.locNy)];
	}

	if ((j == 0) && ((def.source) > 0))
	{
		HostArraysPtr.S_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = def.S_w_gr;
		HostArraysPtr.S_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = def.S_n_gr;
	}

}

void Border_P(ptr_Arrays HostArraysPtr, int i, int j, int k, consts def)
{
	int i1 = i, j1 = j, k1 = k;

	if (i == 0)
	{
		i1 ++;
	}
	if (i == (def.locNx) - 1)
	{
		i1 --;
	}
	if (j == 0)
	{
		j1 ++;
	}
	if (j == (def.locNy) - 1)
	{
		j1 --;
	}
	if ((k == 0) && ((def.locNz) > 2))
	{
		k1 ++;
	}
	if ((k == (def.locNz) - 1) && ((def.locNz) > 2))
	{
		k1 --;
	}
/*
// Если отдельно задаем значения на границах через градиент
	if ((j != 0) && (j != (def.locNy) - 1))
	{
		HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = HostArraysPtr.P_w[i1 + j1 * (def.locNx) + k1 * (def.locNx) * (def.locNy)];
	}
	else if (j == 0)
	{
		HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = def.P_atm;
	}
	else
	{
		HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = HostArraysPtr.P_w[i1 + j1 * (def.locNx) + k1 * (def.locNx) * (def.locNy)];// + ro_eff_gdy(HostArraysPtr, i1, j1, k1, def);
	}
*/
// Если отдельно задаем фиксированные значения на границах

	if ((j != 0) && (j != (def.locNy) - 1))
	{
		HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = HostArraysPtr.P_w[i1 + j1 * (def.locNx) + k1 * (def.locNx) * (def.locNy)];
	}
	else if (j == 0)
	{
		HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = def.P_atm;
	}
	else
	{
		HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = 1.5 * def.P_atm;
	}
/*
// Если внешний слой дублирует первый внутренний
	HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = HostArraysPtr.P_w[i1 + j1 * (def.locNx) + k1 * (def.locNx) * (def.locNy)];
*/
}

void data_initialization(ptr_Arrays HostArraysPtr, long int* t, consts def)
{
	*t = 0;
	for (int i = 0; i < def.locNx; i++)
		for (int j = 0; j < def.locNy; j++)
			for (int k = 0; k < def.locNz; k++)
				if (is_active_point(i, j, k, def))
				{
					// Преобразование локальных координат процессора к глобальным
					int I = local_to_global(i, 'x', def);

					int media = HostArraysPtr.media[i + j * def.locNx + k * def.locNx * def.locNy] = 0;

					// Линейное изменение насыщенностей в начальном распределении
	/*				int j1 = def.locNy / 2;

					if (j < j1)
					{
						HostArraysPtr.S_w[i + j * def.locNx + k * def.locNx * def.locNy] = def.S_w_gr + (def.S_w_init - def.S_w_gr) * j / j1;
						HostArraysPtr.S_n[i + j * def.locNx + k * def.locNx * def.locNy] = def.S_n_gr + (def.S_n_init - def.S_n_gr) * j / j1;
					}
					else
					*/
					if ((j == 0) && ((def.source) > 0))
					{
						HostArraysPtr.S_w[i + j * def.locNx + k * def.locNx * def.locNy] = def.S_w_gr;
						HostArraysPtr.S_n[i + j * def.locNx + k * def.locNx * def.locNy] = def.S_n_gr;
					}
					else
					{
						HostArraysPtr.S_w[i + j * def.locNx + k * def.locNx * def.locNy] = def.S_w_init;
						HostArraysPtr.S_n[i + j * def.locNx + k * def.locNx * def.locNy] = def.S_n_init;
					}

					// Если отдельно задаем значения на границах через градиент
/*					if (j == 0)
					{
						HostArraysPtr.P_w[i + j * def.locNx + k * def.locNx * def.locNy] = def.P_atm;
					}
					else
					{
						HostArraysPtr.P_w[i + j * def.locNx + k * def.locNx * def.locNy] = HostArraysPtr.P_w[i + (j - 1) * def.locNx + k * def.locNx * def.locNy] + ro_eff_gdy(HostArraysPtr, i, j - 1, k, def);
					}
*/
					// Если отдельно задаем фиксированные значения на границах

					if ((j != 0) && (j != (def.locNy) - 1))
					{
						HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = HostArraysPtr.P_w[i + (j - 1) * def.locNx + k * def.locNx * def.locNy] + ro_eff_gdy(HostArraysPtr, i, j - 1, k, def);
					}
					else if (j == 0)
					{
						HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = def.P_atm;
					}
					else
					{
						HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = 1.5 * def.P_atm;
					}
/*
					// Если внешний слой дублирует первый внутренний
					if (j <= 1)
					{
						HostArraysPtr.P_w[i + j * def.locNx + k * def.locNx * def.locNy] = def.P_atm;
					}
					else if (j < (def.Ny) - 1)
					{
						HostArraysPtr.P_w[i + j * def.locNx + k * def.locNx * def.locNy] = HostArraysPtr.P_w[i + (j - 1) * def.locNx + k * def.locNx * def.locNy] + ro_eff_gdy(HostArraysPtr, i, j - 1, k, def);
					}
					else
					{
						HostArraysPtr.P_w[i + j * def.locNx + k * def.locNx * def.locNy] = HostArraysPtr.P_w[i + (j - 1) * def.locNx + k * def.locNx * def.locNy];
					}
*/
					HostArraysPtr.ro_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = def.ro0_w * (1. + (def.beta_w) * (HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] - def.P_atm));

					///!!!! Не учитываются каппилярные силы! Или надо считать перед этим шагом P_w, P_g
					HostArraysPtr.ro_n[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = def.ro0_n * (1. + (def.beta_n) * (HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] - def.P_atm));
					HostArraysPtr.ro_g[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] = def.ro0_g * (1. + (def.beta_g) * (HostArraysPtr.P_w[i + j * (def.locNx) + k * (def.locNx) * (def.locNy)] - def.P_atm));

					test_nan(HostArraysPtr.S_n[i + j * def.locNx + k * def.locNx * def.locNy], __FILE__, __LINE__);
					test_nan(HostArraysPtr.P_w[i + j * def.locNx + k * def.locNx * def.locNy], __FILE__, __LINE__);
					test_nan(HostArraysPtr.media[i + j * def.locNx + k * def.locNx * def.locNy], __FILE__, __LINE__);
					test_nan(HostArraysPtr.S_w[i + j * def.locNx + k * def.locNx * def.locNy], __FILE__, __LINE__);
				}
}

