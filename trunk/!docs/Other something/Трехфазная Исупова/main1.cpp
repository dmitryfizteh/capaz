


//����������� ���������:
#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <stdlib.h>






//-------------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------------
													//���������� ���������:
		//�����:

//����� ������:
#define Version 3001
//������ ���� ������� (�� �����������) � ������:
#define Width 1.8 //0.9	
//������ ���� ������� (�� ���������) � ������:
#define Height 0.64	
//��� �� ������� � ��������:
#define tau 0.005	

//���������� ����� �� ����������� (�� ��� x) (������ ���� �������� ��� ��������� ���������):
#define N1 51 //25	//46	
//���������� ����� �� ��������� (�� ��� y):
#define N2 18	//33	
//���������� ����� �� �������:
#define M 12000000	
//���������� �����, � ������� ����� �������� (��������)
//3
#define N_ist 3

//����� ������� ����� �� ������� ������ ������ � �����:
//1000
#define M_out 10000	

//������������ �� ������������ ������� (1-��, 0-���):
//1
#define Interf_usl 1

//������� ����� (1-����� ����, 0-���):
#define Linz_exist 1
//������� �������� (1-����, 0-���):
#define Zaglushka 1
//������������ �� ��������� ������� S=S0 ������ ��������� (1-��, 0-���):
//1
#define Gr_usl 1
//���� �� ������ �����:
//0
#define Rain 0


//������ ������ ����� �� ���� y � x � ������:
#define lev_y1 0.1875	
#define lev_y2 0.325
#define lev_x1	0.75	//0.1875	
#define lev_x2	1.05	//0.7125


//****************************************************************************
//����� ��������� ���������� ������:

		//��������� ��� ������ �������:
//���������� �������� (������ �� ������������):
#define M_zeid 500	
//�������� �������, ������������ �������� ���������� ������� � ���� ������
#define epsilon 1e-4

		//��������� ������:
//������ 2 ������ ������������ � � ����������
const double m[]={0.4, 0.39};
const double k[]={6.64e-11, 7.15e-15};		// [m^2]
const double S_wr[]={0.09, 0.12};	// eto S_wr
const double lambda[]={2.7, 2.0};
const double p_d[]={755, 2060};		// [Pa]

		//��������� ���������:
const double mu1=0.001;			// [s*Pa]
const double mu2=0.0009;				// [s*Pa]
const double ro1=1000;				// [kg/m^3]
const double ro2=1460;				// [kg/m^3]

		//������ ��������:
//��������� ���������� �������:
const double g=9.8;					// [m/s^2]
//�������� ��������� � ������, ��� �� ����:
//5e-5	0.0
const double Q=5e-5;//5e-5;				// [kg/s]
//���������� �������� S2 �� ������� �������, ���� ������������ ������� ������� ����:
//0.8
const double S0=0.8;				
//����������� ��������:
const double p_atm=0;//100000;			// [Pa]
//****************************************************************************

//� ������ ��� ���������:

//��������� �������� ������������� ���� � �����:
const double S_w_0=0.3;				
const double S_n_0=0.3;				
//���������� �������� �������� ���� (���� ���������):
const double P_g=1e5;
//��������� � ������ �������:
const double alpha = 4.8e-4;
const double sigma_gw = 20.0;
const double sigma_nw = 30.0;
const double sigma_gn = 10.0;
const double n_van_gen = 3.25;
//����������� ���������:
double beta_gn, beta_nw;
//���������� ������������:
const double S_w_r = 0.02;
const double S_g_r = 0.01;
const double S_n_r = 0.01;
		//��������� ���������:
const double mu_w=0.001;			// [s*Pa]
const double mu_n=0.00221;			// [s*Pa]
const double mu_g=0.0000184;		// [s*Pa]
const double ro_w=1000.0;				// [kg/m^3]
const double ro_n=850.0;				// [kg/m^3]
const double ro_g=1.4;				// [kg/m^3]
//���������� �������� S_n �� ������� �������, ���� ������������ ������� ������� ����:
//0.6
const double S_n_ist=0.6;				




		//����������� ���������:
//��� �� �����������:
const double h1=Width/(N1-1);
//��� �� ���������:
const double h2=Height/(N2-1);;
//�������� S �� ������� (����������� � ������� init):
double S_zv;
//������ ������������ ����� �� ���� (����������� � ������� init):
int N_x_1, N_x_2, N_y_1, N_y_2;
//���������� (� ������) ����� ��������� ����� (����������� � ������� init):
double x[N1], y[N2];	
//��������� ��� ��������� ������ ������������ �����:
struct Interf_points
{
	int i,j;
	double x1,y1;
};
//���������� ����� ���� (������������ ���������� ������).
Interf_points *i_p;
//���-�� ��������� � ���� �������:
int N_i_p;



//-------------------------------------------------------------------------------------------------------------------------------
									//����������:
//����� ���� �� �������
int vrema;	


		//�������� ��������:
//��������� ��� ������������ (���������� ����� ���� �������� ������� �� ������ ���� �� �������):
struct Sat	
{
	//�������� ������������ �� �������� (������ 0) � ����������� (������ 1) �����:
	double value[2];
	//������������ ����, ����� � ���� ��������������:
	double w[2];
	double n[2];
	double g[2];
};

//��������� ��� �������� (���������� ����� ���� �������� ������� �� ������ �������� � ������ �������):
struct Pres
{
	//�������� ��������:
//bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb	
	double value;
//eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee
	//�������� ���� � ����� (� ���� - ���������):
	//double w,n;
	double w[2];
	double n[2];
};

//��������� ��� ����������, � ������� ����� ���� ����. � �������� �������� (� ��������, �-�� ������������):
struct Func
{
	//������� �� ������������ � ���������� ������:
	double k1, k2, S_e, p_k, a_zeid;

	//������� �� ������������ � ���������� ������:
	double p_cgn, p_cnw;
	double S_line_w, S_line_n, S_line_g;
	double k_rn, k_rw, k_rg;

	//����� �����:
	int sand;
};

//��������� ��� �������� ������������� (0 - ������., 1 - ����.) � ������ ������ �� �����
struct Virt_type
{
	int l0,r0,u0,d0;
};

//��������� ��� ������������� � ������ ����������:
struct Relax_koef
{
	double znam;
	//����-�� ��� p[i][j-1], p[i][j+1], p[i-1][j], p[i+1][j] �����-��, � ��������� ����-�:
	double verh, niz, levo, pravo,f;
};

//�������� ��������� ��� ������ � �����:
struct Point
{
	//�������� � ����. �������� ���������� � ��������� Func:
	Func f_value[2];
	//����� ����� ����� ������ ���:
	int cl_x, cl_y;
	//������������� � ������ ������:
	Virt_type vrt;
	//������������ � ������ ����������:
	Relax_koef rel;
};


		//���������� ����� ���� ��������:
//������������ �� ������� � ��������� ����� �� �������, � ����� �������� ����������:
Sat **S, **S_next, **S_buf;

//�������� � ������� � ��������� ��������� �������, � ��������:
Pres **p, **p_next, **p_buf;

//������ � ������:
Point pnt[N1][N2];






//-------------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------------
													//��������� �������:

											//����������� �������:
//��������� ������ ��� ������������ �������:
void mem_alloc();
//������������ ������:
void mem_free();
//������������� ����������:
void init();
//���� ��������� ������� ��� ������������ � ��������:
void nach_usl(int, int);
//������������ ������������� (���������� virt) ��� ���� �����:
void virt_init();
//���������� ������ ����� �����:
void sand_init();
//����������� ���� ����� ����� ���������� ���:
int point_class(int i, int os);
//������������� i_p
void init_i_p();
//���������� � init() ��� ���������� �����:
void init_wo_linza();
//���������� � init() ��� ���������� ����������:
void init_wo_interf();
//������ � ���� ��������� ����������:
void rec_info();


											//������������ �������:
//������ ���������� ���� � ����
void rec();
//������ ������� � ����
void rec2();
//������ ������������ ���� � ������������ ������
void rec_virt(FILE *);
//��� �� ������� (������);
void step_time();		
//���������� ������� �� S � ������ ��������������� ������� (�.�. ���������� f_value)
void step_0();
//���������� �������� (������������ � ���� 0):
void step_0_pressure(int i,int j, int v);
//���������� k (������������ � ���� 0):
void step_0_k(int i,int j, int v);
//������ ����� ���� 0:
void step_0_b();
//������ �-�� divgrad() � ����� ������ ������:
double step_0_divgrad(int i,int j);
//������� ��������� ��� ������������:
void step_1();
//������� ��������� ��� ������������ ����:
void step_1_s_w();
//������� ��������� ��� ������������ �����:
void step_1_s_n();
//������ ������������ ����:
void step_1_s_g();
//��������� ������� ��� ������������:
void step_1_gr();
//��������� ������� ��� ��������:
void step_1_gr_pres();
//���������� ������������ ���� ������������ � ���� 1
void step_1_virt();
//��������, ����� �� ����� �� ����������:
int is_interf(int,int);
//��������, ���� �� �������� � ������ �����:
int is_ist(int, int);
//������� ��������� ��� ��������:
/*void step_2();
//���� �������� � ������ ����������:
void step_2_relax_iter();
//�������, ��������� ������� � ��-�� ��� ��������:
double step_2_nevazka();
//��������� ������� ��� �������� (������� ����������� �� ��������):
void step_2_gr_before();
//��������� ������� ��� �������� (������� ����������� ����� ��������):
void step_2_gr_after();
//��� ������ ��������� ��������:
void step_2_gr_after_2();*/
//���������� ��������� � ����:
void save();
//�������� ��������� �� �����:
void load();




											//��������������� �������:
//��������������� �-�� ��� �������������� �������� ����������:
void test_print(int); 
//��������������� �-�� ��� ������ � ����� ����� ������: 
void test_rec(int);
//�-�� ��� ������������ ������� init
void test1();
//��� ������������ � ������� ����������� �������:
void test_step_1();
//abs ��� ����� double
double absd(double);
											//���������������� �-��:
//�-� ��� ������ � ����� � 2-������, 3-������ ����������� ���� (������ �� ������������):
void bu_rec();		



//-------------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------------
													//�������:


main()
{
	mem_alloc();
	init();
	rec_info();
	//test_print(103);
	//rec();
	//test_rec(15);
	//test1();
	//printf("%f\n",Q);

	vrema=0;

	//load();

	//test_rec(2);

	for(;vrema<M;vrema++)
	{
		step_time();
	}
	
	mem_free();

	return 0;
}

//-------------------------------------------------------------------------------------------------------------------------------

void mem_alloc()
{
	int i;
	//���� ������:
	int err=0;

	//������ ��� ������ S:
	S=(Sat **)malloc(N1*sizeof(Sat *)); 
    for(i=0;i<N1;i++) 	S[i]=(Sat *)malloc((N2)*sizeof(Sat));		

	S_next=(Sat **)malloc(N1*sizeof(Sat *)); 
    for(i=0;i<N1;i++) S_next[i]=(Sat *)malloc((N2)*sizeof(Sat));	


	//������ ��� ������ p:
	p=(Pres **)malloc(N1*sizeof(Pres *)); 
    for(i=0;i<N1;i++) p[i]=(Pres *)malloc((N2)*sizeof(Pres));

	p_next=(Pres **)malloc(N1*sizeof(Pres *)); 
    for(i=0;i<N1;i++) p_next[i]=(Pres *)malloc((N2)*sizeof(Pres));

	
	//�������� �� ������:
	if(S==NULL || S_next==NULL || p==NULL || p_next==NULL) err=1;
	for(i=0;i<N1;i++)
	{
		if(S[i]==NULL || S_next[i]==NULL || p[i]==NULL || p_next[i]==NULL) err=1;
	}

	if(err==1)
	{
		printf("������ ��������� ������");
		getchar();
	}

}

//-------------------------------------------------------------------------------------------------------------------------------

void mem_free()
{
	int i;

	
	for(i=0;i<N1;i++)
	{
		free(S[i]);
		free(S_next[i]);
		free(p[i]);
		free(p_next[i]);
	}

	free(S);
	free(S_next);
	free(p);
	free(p_next);
	
}

//-------------------------------------------------------------------------------------------------------------------------------

void init()
{
	int i, j;

//bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb
	/*
	//���������� ����������� ��������� S �� �������:
	S_zv=(1-S_wr[0])*powf(p_d[0]/p_d[1],lambda[0])+S_wr[0];
	printf("S_zv = %f\n",S_zv);
	*/
//eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee	

	vrema=0;

	//������ ��� ��� ������ �������:
	beta_gn = sigma_gw / sigma_gn;
	beta_nw = sigma_gw / sigma_nw;

	//������ ��������� �� ������ �����:
	for(i=0;i<N1;i++)
		x[i]=i*h1;
	for(j=0;j<N2;j++)
		y[j]=j*h2;

	//������ ������� N_x_1, N_x_2, N_y_1, N_y_2 (����� ��������� � ��������, ���� �� �������� ����� � �����, �.�. ���������):
	N_x_1=lev_x1/h1;
	N_x_2=lev_x2/h1;
	N_y_1=lev_y1/h2;
	N_y_2=lev_y2/h2;
	if(x[N_x_2]<lev_x2) N_x_2++;
	if(y[N_y_2]<lev_y2) N_y_2++;
	//�������� �� ������������ ������:
	if(N_x_1<=0 || N_y_1<=0 || N_x_2<=0 || N_y_2<=0 || N_x_1>=N1 || N_y_1>=N2 || N_x_2>=N1 || N_y_2>=N2)
	{
		printf("������. ������ ������ ����� �� �������� � ��������� �������\n");
		getchar();
	}
	


	//������������� �������, ��������� �� �����
	for(i=0;i<N1;i++)
	{
		for(j=0;j<N2;j++)
		{
			//��������� ������� ��� ������������ � ��������:
			nach_usl(i,j);
			pnt[i][j].cl_x=point_class(i,1);
			pnt[i][j].cl_y=point_class(j,2);
		}
	}


	virt_init();
	sand_init();
	//�� ������������:
	//init_i_p();

	if(Linz_exist==0) init_wo_linza();
	if(Interf_usl==0) init_wo_interf();
}

//-------------------------------------------------------------------------------------------------------------------------------

void nach_usl(int i, int j)
{
	//������� ������� ��� ������������:
//bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb
	/*
	S[i][j].value[0]=0;
	S[i][j].value[1]=0;
	S_next[i][j].value[0]=0;
	S_next[i][j].value[1]=0;
	*/
//eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee

	S[i][j].w[0]=S_w_0;
	S[i][j].w[1]=S_w_0;
	S[i][j].n[0]=S_n_0;
	S[i][j].n[1]=S_n_0;
	S[i][j].g[0]=1-S_w_0-S_n_0;
	S[i][j].g[1]=1-S_w_0-S_n_0;

	S_next[i][j].w[0]=S_w_0;
	S_next[i][j].w[1]=S_w_0;
	S_next[i][j].n[0]=S_n_0;
	S_next[i][j].n[1]=S_n_0;
	S_next[i][j].g[0]=1-S_w_0-S_n_0;
	S_next[i][j].g[1]=1-S_w_0-S_n_0;

//bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb
	/*
	//�������� - ����������������:
	p[i][j].value=p_atm+ro1*g*j*h2;
	p_next[i][j].value=p[i][j].value;


	//���������� �������� ������ � ��������� ������ (��� �������� � ������ ����������) ��� �������,
	//� ������� � ��������� ������ �������� �� ������������ (���� �������� ���� ��� ����� ������):
		pnt[i][j].rel.f=0;
		pnt[i][j].rel.levo=0;
		pnt[i][j].rel.niz=0;
		pnt[i][j].rel.pravo=0;
		pnt[i][j].rel.verh=0;
		pnt[i][j].rel.znam=0;
*/
//eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee
		
}

//-------------------------------------------------------------------------------------------------------------------------------

//���� ����� (��� �������� ��������� ��� ��� x, � ��� y - ����������):
//0, 1 - ����� � ������ ������� ������� �����-��
//2, 3 - ���������� ����� �������� ������� � ������ �����-�� ����� � ������
//4, 5 - �����-�� ����� � ������ ���������� �����
//6 - ���������� ������� �����
//�.�. ������� �����: 0246531
int point_class(int i, int os)
{
	//i - ����� ����� ����� ��������������� ���; os - ����� ��� (���� 1, �� x, ���� 2, �� y)

	//������ ������ � ������ ������� ������� ������
	int gr1, gr2;
	//���-�� ����� ����� �����. ���
	int N_max;

	//����������� ���������� � ����������� �� ���:
	if(os==1)
	{
		gr1=N_x_1;
		gr2=N_x_2;
		N_max=N1;
	}
	else if(os==2)
	{
		gr1=N_y_1;
		gr2=N_y_2;
		N_max=N2;
	}
	else 
	{
		printf("������ ������������� ������� point_class. ����������� ������� ���.");
		getchar();
	}


	//����������� ���� �����:
	if(i==0) return 0;
	else if(i==N_max-1) return 1;
	else if(i>0 && i<gr1) return 2;
	else if(i<N_max-1 && i>gr2) return 3;
	else if(i==gr1) return 4;
	else if(i==gr2) return 5;
	else if(i>gr1 && i<gr2) return 6;
	else 
	{
		printf("������ ������������� ������� point_class. ����� ������� �� ������� �������.");
		getchar();
	}


	return 10;
}

//-------------------------------------------------------------------------------------------------------------------------------

void virt_init()
{
	int i, j;

	//������� �� �������� ������, � �����, ��� ����, �������� �������:
	for(i=0;i<N1;i++)
	{
		for(j=0;j<N2;j++)
		{
			pnt[i][j].vrt.r0=0;
			pnt[i][j].vrt.l0=0;
			pnt[i][j].vrt.u0=0;
			pnt[i][j].vrt.d0=0;
		}
	}

	for(i=0;i<N1;i++)
	{
		for(j=0;j<N2;j++)
		{
//*********
	if(pnt[i][j].cl_x>1 && pnt[i][j].cl_y>1)
	{
		//����� ��� x:
		//����� ������� �����:
		if(pnt[i][j].cl_y==6 && pnt[i][j].cl_x==4)
		{
			pnt[i][j].vrt.r0=1;
			pnt[i+1][j].vrt.l0=1;
		}
		//������ �������:
		if(pnt[i][j].cl_y==6 && pnt[i][j].cl_x==5)
		{
			pnt[i][j].vrt.l0=1;
			pnt[i-1][j].vrt.r0=1;
		}
		//����� ��� y:
		//������� ������� �����:
		if(pnt[i][j].cl_x==6 && pnt[i][j].cl_y==4)
		{
			pnt[i][j].vrt.d0=1;
			pnt[i][j+1].vrt.u0=1;
		}
		//������ �������:
		if(pnt[i][j].cl_x==6 && pnt[i][j].cl_y==5)
		{
			pnt[i][j].vrt.u0=1;
			pnt[i][j-1].vrt.d0=1;
		}
	}
//*********
		}
	}

}

//-------------------------------------------------------------------------------------------------------------------------------

void sand_init()
{
	int i, j;

	for(i=0;i<N1;i++)
	{
		for(j=0;j<N2;j++)
		{
			//���� ��� �����:
			{
				pnt[i][j].f_value[0].sand=0;
				pnt[i][j].f_value[1].sand=0;
			}

			//���� ������ �����:
			if(pnt[i][j].cl_x==6 && pnt[i][j].cl_y==6)
			{
				pnt[i][j].f_value[0].sand=1;				
			}
			//���� �� ���������� ��� ������:
			if(pnt[i][j].cl_x>=4 && pnt[i][j].cl_y>=4)
			{				
				pnt[i][j].f_value[1].sand=1;
			}
						
		}
	}
}

//-------------------------------------------------------------------------------------------------------------------------------

void init_wo_linza()
{
	int i, j;

	//������������� �������, ��������� �� �����
	for(i=0;i<N1;i++)
	{
		for(j=0;j<N2;j++)
		{	
			//���������. ����� ����
			if(pnt[i][j].cl_x>1) pnt[i][j].cl_x=2;
			if(pnt[i][j].cl_y>1) pnt[i][j].cl_y=2;

			//�������������:
			pnt[i][j].vrt.r0=0;
			pnt[i][j].vrt.l0=0;
			pnt[i][j].vrt.u0=0;
			pnt[i][j].vrt.d0=0;

			//����� �����:
			{
				pnt[i][j].f_value[0].sand=0;
				pnt[i][j].f_value[1].sand=0;
			}

		}
	}

	
}

//-------------------------------------------------------------------------------------------------------------------------------

void init_wo_interf()
{
	int i, j;

	//������������� �������, ��������� �� �����
	for(i=0;i<N1;i++)
	{
		for(j=0;j<N2;j++)
		{
			//�������������:
			pnt[i][j].vrt.r0=0;
			pnt[i][j].vrt.l0=0;
			pnt[i][j].vrt.u0=0;
			pnt[i][j].vrt.d0=0;
		
		}
	}

}

//-------------------------------------------------------------------------------------------------------------------------------

void init_i_p()
{
	int i,j, numb;
	const double dobavka_x=h1/10;
	const double dobavka_y=h2/10;

	//���������� N_i_p:
	numb=0;
	for(i=0;i<N1;i++)
	{
		for(j=0;j<N2;j++)
		{
			if(is_interf(i,j))
			{
				numb++;			
			}
		}
	}
	N_i_p=numb;

	//��������� ������ ��� �������:
	i_p=(Interf_points *)malloc((N_i_p)*sizeof(Interf_points));		

	//���������� �������:
	numb=0;
	for(i=0;i<N1;i++)
	{
		for(j=0;j<N2;j++)
		{
			if(is_interf(i,j))
			{
				i_p[numb].i=i;
				i_p[numb].j=j;

					//���������� x1:
				//����� �������:
				if(pnt[i][j].cl_x==4) i_p[numb].x1=x[i]+dobavka_x;
				//������� �����:
				else if(pnt[i][j].cl_x==6) i_p[numb].x1=x[i];
				//������ �������:
				else if(pnt[i][j].cl_x==5) i_p[numb].x1=x[i]-dobavka_x;
				else 
				{
					printf("������ � ������� init_i_p()\n");
					getchar();
				}

					//���������� y1:
				//������� �������:
				if(pnt[i][j].cl_y==4) i_p[numb].y1=y[j]+dobavka_y;
				//������� �����:
				else if(pnt[i][j].cl_y==6) i_p[numb].y1=y[j];
				//������ �������:
				else if(pnt[i][j].cl_y==5) i_p[numb].y1=y[j]-dobavka_y;
				else 
				{
					printf("������ � ������� init_i_p()\n");
					getchar();
				}

				numb++;
			}
		}
	}
}

//-------------------------------------------------------------------------------------------------------------------------------

void rec()
{
	//��������� �� �����:
	FILE *f3S, *f3p;

	FILE *fS_w, *fS_n, *fS_g;
	FILE *fp_w, *fp_n;

	int numb,j,i;
	//������ ��� ����� �����:
	char namefl[30];
	
	numb=vrema/M_out;

					//��������� ����� ������ � ������� ��, � ��������� � �����������:
	sprintf(namefl, "output\\S_w\\Sw_%d.dat",vrema+10000000);
    fS_w=fopen(namefl,"w");
	sprintf(namefl, "output\\S_n\\Sn_%d.dat",vrema+10000000);
    fS_n=fopen(namefl,"w");
	sprintf(namefl, "output\\S_g\\Sg_%d.dat",vrema+10000000);
    fS_g=fopen(namefl,"w");

	sprintf(namefl, "output\\p_w\\Sw_%d.dat",vrema+10000000);
    fp_w=fopen(namefl,"w");
	sprintf(namefl, "output\\p_n\\Sn_%d.dat",vrema+10000000);
    fp_n=fopen(namefl,"w");


//bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb
	/*
						//��������� ����� ������ � ������� ��, � ��������� � �����������:
	sprintf(namefl, "out1s\\S2_%d.dat",vrema);
    f3S=fopen(namefl,"w");
	sprintf(namefl, "out1p\\p1_%d.dat",vrema);
    f3p=fopen(namefl,"w");


					//������ � �����:	
	//������ �����:
	fprintf(f3S,"TITLE =  \"Saturation\" \n");
	fprintf(f3S,"VARIABLES = \"X\",\"Y\",\"S2\" \n");
    fprintf(f3S,"ZONE T = \"BIG ZONE\", J=%d  ,I=%d, F = POINT \n", N2, N1);

	fprintf(f3p,"TITLE =  \"Pressure\" \n");
	fprintf(f3p,"VARIABLES = \"X\",\"Y\",\"P1\" \n");
    fprintf(f3p,"ZONE T = \"BIG ZONE\", J=%d  ,I=%d, F = POINT \n", N2, N1);

	//������ ������:
	for(j=0; j<N2; j++)
	{
			for(i=0; i<N1; i++)
			{
			     fprintf(f3S,"%6.4e %6.4e %6.4e \n",x[i],y[j],S[i][j].value[0]);
				 fprintf(f3p,"%6.4e %6.4e %6.4e \n",x[i],y[j],p[i][j].value);
			}
	 }

	//���������� � ���� � ������������� ����������� �������� � ������������ ������
	//rec_virt(f3S);

	fclose(f3S);
	fclose(f3p);
*/
//eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee

					//������ � �����:	
	//������ �����:
	fprintf(fS_w,"TITLE =  \"Saturation\" \n");
	fprintf(fS_w,"VARIABLES = \"X\",\"Y\",\"Sw\" \n");
    fprintf(fS_w,"ZONE T = \"BIG ZONE\", J=%d  ,I=%d, F = POINT \n", N2, N1);

	fprintf(fS_n,"TITLE =  \"Saturation\" \n");
	fprintf(fS_n,"VARIABLES = \"X\",\"Y\",\"Sn\" \n");
    fprintf(fS_n,"ZONE T = \"BIG ZONE\", J=%d  ,I=%d, F = POINT \n", N2, N1);

	fprintf(fS_g,"TITLE =  \"Saturation\" \n");
	fprintf(fS_g,"VARIABLES = \"X\",\"Y\",\"Sg\" \n");
    fprintf(fS_g,"ZONE T = \"BIG ZONE\", J=%d  ,I=%d, F = POINT \n", N2, N1);


	fprintf(fp_w,"TITLE =  \"Pressure\" \n");
	fprintf(fp_w,"VARIABLES = \"X\",\"Y\",\"Pw\" \n");
    fprintf(fp_w,"ZONE T = \"BIG ZONE\", J=%d  ,I=%d, F = POINT \n", N2, N1);

	fprintf(fp_n,"TITLE =  \"Pressure\" \n");
	fprintf(fp_n,"VARIABLES = \"X\",\"Y\",\"Pn\" \n");
    fprintf(fp_n,"ZONE T = \"BIG ZONE\", J=%d  ,I=%d, F = POINT \n", N2, N1);


	//������ ������:
	for(j=0; j<N2; j++)
	{
			for(i=0; i<N1; i++)
			{
			     fprintf(fS_w,"%6.4e %6.4e %6.4e \n",x[i],y[j],S[i][j].w[0]);
				 fprintf(fS_n,"%6.4e %6.4e %6.4e \n",x[i],y[j],S[i][j].n[0]);
				 fprintf(fS_g,"%6.4e %6.4e %6.4e \n",x[i],y[j],S[i][j].g[0]);

				 fprintf(fp_w,"%6.4e %6.4e %6.4e \n",x[i],y[j],p[i][j].w[0]);
				 fprintf(fp_n,"%6.4e %6.4e %6.4e \n",x[i],y[j],p[i][j].n[0]);
			}
	 }



	
	//�������� ���� �������� ������:
	fclose(fS_w);
	fclose(fS_n);
	fclose(fS_g);

	fclose(fp_w);
	fclose(fp_n);

}

//-------------------------------------------------------------------------------------------------------------------------------

void rec2()
{
	FILE *fS_w, *fS_n, *fS_g;
	//FILE *fp_w, *fp_n;

	int numb,j,i;
	//������ ��� ����� �����:
	char namefl[30];
	
	numb=vrema/M_out;

					//��������� ����� ������ � ������� ��, � ��������� � �����������:
	sprintf(namefl, "output2\\S_w\\Sw_%d.dat",vrema+10000000);
    fS_w=fopen(namefl,"w");
	sprintf(namefl, "output2\\S_n\\Sn_%d.dat",vrema+10000000);
    fS_n=fopen(namefl,"w");
	sprintf(namefl, "output2\\S_g\\Sg_%d.dat",vrema+10000000);
    fS_g=fopen(namefl,"w");

	/*sprintf(namefl, "output2\\p_w\\Sw_%d.dat",vrema+10000000);
    fp_w=fopen(namefl,"w");
	sprintf(namefl, "output2\\p_n\\Sn_%d.dat",vrema+10000000);
    fp_n=fopen(namefl,"w");*/


						//������ � �����:

	for(j=0;j<N2;j++)
	{
		//l=sand(j);
		fprintf(fS_w,"%6.4e ",y[j]);
		fprintf(fS_w,"%8.6e\n",S[N2/2][j].w[0]);
		fprintf(fS_n,"%6.4e ",y[j]);
		fprintf(fS_n,"%8.6e\n",S[N2/2][j].n[0]);
		fprintf(fS_g,"%6.4e ",y[j]);
		fprintf(fS_g,"%8.6e\n",S[N2/2][j].g[0]);
		//fprintf(fp2,"%6.4e ",y[j]);
		//fprintf(fp2,"%8.6e\n",p1_a[1][j]+re[1][j].p_k);
	}



	//�������� ���� �������� ������:
	fclose(fS_w);
	fclose(fS_n);
	fclose(fS_g);


}

//-------------------------------------------------------------------------------------------------------------------------------

void rec_virt(FILE *fpt)
{
	int ii;
	
	//������ ������:
	for(ii=0; ii<N_i_p; ii++)
	{
			//fprintf(fpt,"%6.4e %6.4e %6.4e \n",i_p[ii].x1,i_p[ii].y1,S[i_p[ii].i][i_p[ii].j].value[1]);
			fprintf(fpt,"%6.4e %6.4e %6.4e \n",i_p[ii].x1,i_p[ii].y1,1.0);
	}
	
}

//-------------------------------------------------------------------------------------------------------------------------------

void step_time()
{
	step_0();
	//� ����� ������ ���� ��������� �������� ���������, ����� ��������������� ���. ���-� ������������:
	//if(vrema==0) step_2();
	step_1();
	//test_step_1();
	//step_2();

	//������ ������� ���� ������������:
	S_buf=S;
	S=S_next;
	S_next=S_buf;


	//����� �� ����� ������ ���� �� �������:
	printf("vrema = %i;\n",vrema);

	//������ � ����� ��� ���������� � ���������� ���������:
	if(vrema%M_out==0)
	{
		printf("attantion, I'm writing, vrema = %d\n",vrema);
		rec();
		rec2();
		save();
		printf("ok\n");
	}
}

//-------------------------------------------------------------------------------------------------------------------------------

void step_0()
{
	//����������, ����� ����� � �������������
	int i, j, sn, v;
	//��������� �������� S_e (��� ���������� ������� ������):
	double loc_S_e;


	//������� ��������, ��������� �� �������� � ����� �����:
	for(i=0;i<N1;i++)
	{
		for(j=0;j<N2;j++)
		{
			for(v=0;v<2;v++)
			{
				//����� ����� � ������ �����:
				sn=pnt[i][j].f_value[v].sand;

				step_0_pressure(i, j, v);
				step_0_k(i, j, v);

//bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb
				/*
				//�������� S_e:
				loc_S_e=(1-S[i][j].value[v]-S_wr[sn])/(1-S_wr[sn]);
				pnt[i][j].f_value[v].S_e=loc_S_e;
				//�������� k1:
				pnt[i][j].f_value[v].k1=powf(loc_S_e,(2+3*lambda[sn])/lambda[sn]);
				//�������� k2:
				pnt[i][j].f_value[v].k2=powf(1-loc_S_e,2)*(1-powf(loc_S_e,(2+lambda[sn])/lambda[sn]));
				//�������� p_k:
				pnt[i][j].f_value[v].p_k=p_d[sn]*powf(loc_S_e,-(1/lambda[sn]));
				//�������� a_zeid:
				pnt[i][j].f_value[v].a_zeid=k[sn]*(pnt[i][j].f_value[v].k1/mu1+pnt[i][j].f_value[v].k2/mu2);
				*/
//eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee
			}
		}
	}


	step_1_gr_pres();

//bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb
	//������ ��������, ��������� �� ���������� �����:
	//step_0_b();
//eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee

}

//-------------------------------------------------------------------------------------------------------------------------------

void step_0_b()
{
	int i, j, sn, sn_u, sn_d;

	//��������� �������������:
	struct 
	{
		int l, r, u, d;
	} v;

	//��������� ����� ��� f:
	double vsp1, vsp2, vsp3;

	for(i=1;i<N1-1;i++)
	{
		for(j=1;j<N2-1;j++)
		{
			//���������� v:
			v.l=pnt[i][j].vrt.l0;
			v.r=pnt[i][j].vrt.r0;
			v.u=pnt[i][j].vrt.u0;
			v.d=pnt[i][j].vrt.d0;

			//����. ��� p[i][j]:
			pnt[i][j].rel.znam=(pnt[i+1][j].f_value[v.r].a_zeid+pnt[i][j].f_value[v.r].a_zeid+pnt[i][j].f_value[v.l].a_zeid+pnt[i-1][j].f_value[v.l].a_zeid)/(2*h1*h1)+(pnt[i][j+1].f_value[v.d].a_zeid+pnt[i][j].f_value[v.d].a_zeid+pnt[i][j].f_value[v.u].a_zeid+pnt[i][j-1].f_value[v.u].a_zeid)/(2*h2*h2);
			//��������� ����-��:
			pnt[i][j].rel.pravo=(pnt[i+1][j].f_value[v.r].a_zeid+pnt[i][j].f_value[v.r].a_zeid)/(2*h1*h1);
			pnt[i][j].rel.levo=(pnt[i][j].f_value[v.l].a_zeid+pnt[i-1][j].f_value[v.l].a_zeid)/(2*h1*h1);
			pnt[i][j].rel.niz=(pnt[i][j+1].f_value[v.d].a_zeid+pnt[i][j].f_value[v.d].a_zeid)/(2*h2*h2);
			pnt[i][j].rel.verh=(pnt[i][j].f_value[v.u].a_zeid+pnt[i][j-1].f_value[v.u].a_zeid)/(2*h2*h2);

			//������ ������� f:
			sn=pnt[i][j].f_value[v.u].sand;
			vsp1=ro1*g/mu1*(pnt[i][j].f_value[v.u].k1* k[sn]-pnt[i][j-1].f_value[v.u].k1*k[sn])/h2;
			vsp2=ro2*g/mu2*(pnt[i][j].f_value[v.u].k2* k[sn]-pnt[i][j-1].f_value[v.u].k2*k[sn])/h2;
			vsp3=step_0_divgrad(i,j)/mu2;


			//temp
			/*sn_d=pnt[i][j+1].f_value[v.d].sand;
			sn_u=pnt[i][j-1].f_value[v.u].sand;
			vsp1=ro1*g/mu1*(pnt[i][j+1].f_value[v.d].k1* k[sn_d]-pnt[i][j-1].f_value[v.u].k1*k[sn_u])/(2*h2);
			vsp2=ro2*g/mu2*(pnt[i][j+1].f_value[v.d].k2* k[sn_d]-pnt[i][j-1].f_value[v.u].k2*k[sn_u])/(2*h2);
			//temp

			//temp
			//v.d=pnt[i][j+1].vrt.u0;
			v.d=0;
			v.u=0;
			sn_d=pnt[i][j].f_value[v.d].sand;
			sn_u=pnt[i][j-1].f_value[v.u].sand;
			vsp1=ro1*g/mu1*(pnt[i][j].f_value[v.d].k1* k[sn_d]-pnt[i][j-1].f_value[v.u].k1*k[sn_u])/(h2);
			vsp2=ro2*g/mu2*(pnt[i][j].f_value[v.d].k2* k[sn_d]-pnt[i][j-1].f_value[v.u].k2*k[sn_u])/(h2);*/
			//temp

			/*if(Interf_usl==0)
			{
				sn=pnt[i][j].f_value[v.u].sand;
				sn_u=pnt[i][j-1].f_value[v.u].sand;
				vsp1=ro1*g/mu1*(pnt[i][j].f_value[v.u].k1* k[sn]-pnt[i][j-1].f_value[v.u].k1*k[sn_u])/h2;
				vsp2=ro2*g/mu2*(pnt[i][j].f_value[v.u].k2* k[sn]-pnt[i][j-1].f_value[v.u].k2*k[sn_u])/h2;
			}*/

			pnt[i][j].rel.f=-vsp1-vsp2+vsp3;
			//���� ���������:
			if(is_ist(i,j)) pnt[i][j].rel.f+=Q/(h2*h2);

		}
	}


}

//-------------------------------------------------------------------------------------------------------------------------------

double step_0_divgrad(int i,int j)
{
	//����������� div(k*k2* grad(P_k))

	double b_x2, b_x1, b_y2, b_y1;
	double ax2, ax1, ay2, ay1;
	
	//��������� �������������:
	struct 
	{
		int l, r, u, d;
	} v;


	//���������� v:
			v.l=pnt[i][j].vrt.l0;
			v.r=pnt[i][j].vrt.r0;
			v.u=pnt[i][j].vrt.u0;
			v.d=pnt[i][j].vrt.d0;


	//l=sand(i,j);

	b_x2=(pnt[i+1][j].f_value[v.r].p_k-pnt[i][j].f_value[v.r].p_k)/h1;
	b_x1=(pnt[i][j].f_value[v.l].p_k-pnt[i-1][j].f_value[v.l].p_k)/h1;
	b_y2=(pnt[i][j+1].f_value[v.d].p_k-pnt[i][j].f_value[v.d].p_k)/h2;
	b_y1=(pnt[i][j].f_value[v.u].p_k-pnt[i][j-1].f_value[v.u].p_k)/h2;
	ax2=(pnt[i+1][j].f_value[v.r].k2+pnt[i][j].f_value[v.r].k2)*k[pnt[i][j].f_value[v.r].sand]/2;
	ax1=(pnt[i][j].f_value[v.l].k2+pnt[i-1][j].f_value[v.l].k2)*k[pnt[i][j].f_value[v.l].sand]/2;
	ay2=(pnt[i][j+1].f_value[v.d].k2+pnt[i][j].f_value[v.d].k2)*k[pnt[i][j].f_value[v.d].sand]/2;
	ay1=(pnt[i][j].f_value[v.u].k2+pnt[i][j-1].f_value[v.u].k2)*k[pnt[i][j].f_value[v.u].sand]/2;
	return (ax2*b_x2-ax1*b_x1)/h1+(ay2*b_y2-ay1*b_y1)/h2;

}

//-------------------------------------------------------------------------------------------------------------------------------

void step_0_pressure(int i,int j, int v)
{
	//��, ��� � ���������� �������:
	double tmp1, tmp2, sum;

	//������ S � ������:
	sum = S_w_r + S_n_r + S_g_r;
	pnt[i][j].f_value[v].S_line_w = (S[i][j].w[v]-S_w_r) / (1-sum);
	pnt[i][j].f_value[v].S_line_n = (S[i][j].n[v]-S_n_r) / (1-sum);
	pnt[i][j].f_value[v].S_line_g = (S[i][j].g[v]-S_g_r) / (1-sum);

	//������ p_cgn � p_cnw:
	tmp1 = powf(pnt[i][j].f_value[v].S_line_w, n_van_gen/(1-n_van_gen)) - 1;
	pnt[i][j].f_value[v].p_cnw = 1/(alpha*beta_nw) * powf(tmp1, 1/n_van_gen);

	tmp2 = powf(1-pnt[i][j].f_value[v].S_line_g, n_van_gen/(1-n_van_gen)) - 1;
	pnt[i][j].f_value[v].p_cgn = 1/(alpha*beta_gn) * powf(tmp2, 1/n_van_gen);


	//������ p_n � p_w:
	p[i][j].n[v]=P_g-pnt[i][j].f_value[v].p_cgn;
	p[i][j].w[v]=P_g-(pnt[i][j].f_value[v].p_cgn+pnt[i][j].f_value[v].p_cnw);
}

//-------------------------------------------------------------------------------------------------------------------------------

//��� ������� ������ ����������� ������ ����� step_0_pressure
void step_0_k(int i,int j, int v)
{
	double  k_rnw, k_rng;
	double tmp_rng1, tmp_rng2, tmp_rnw;
	double tmp_rw1, tmp_rw2, tmp_rg;

	//������ k_rn:
		//������ k_rnw:
		tmp_rnw = 1 - powf(pnt[i][j].f_value[v].S_line_w, n_van_gen/(n_van_gen-1));
		k_rnw = powf(1-pnt[i][j].f_value[v].S_line_w, 0.5) * powf(tmp_rnw, 2*(n_van_gen-1)/n_van_gen);
		//������ k_rng:
		tmp_rng2 = 1 - powf(pnt[i][j].f_value[v].S_line_n, n_van_gen/(n_van_gen-1));
		tmp_rng1 = 1 - powf(tmp_rng2, (n_van_gen-1)/n_van_gen);
		k_rng = powf(pnt[i][j].f_value[v].S_line_n, 0.5) * powf(tmp_rng1, 2);
	pnt[i][j].f_value[v].k_rn = pnt[i][j].f_value[v].S_line_n * k_rnw * k_rng / ((1-S[i][j].w[v])*(S[i][j].w[v]+S[i][j].n[v]));

	//������ k_rw:
	tmp_rw2 = 1 - powf(pnt[i][j].f_value[v].S_line_w, n_van_gen/(n_van_gen-1));
	tmp_rw1 = 1 - powf(tmp_rw2, (n_van_gen-1)/n_van_gen);
	pnt[i][j].f_value[v].k_rw = powf(pnt[i][j].f_value[v].S_line_w, 0.5) * powf(tmp_rw1, 2);

	//������ k_rg:
	tmp_rg = 1 - powf(1-pnt[i][j].f_value[v].S_line_g, n_van_gen/(n_van_gen-1));
	pnt[i][j].f_value[v].k_rg = powf(pnt[i][j].f_value[v].S_line_g, 0.5) * powf(tmp_rg, 2*(n_van_gen-1)/n_van_gen);
}

//-------------------------------------------------------------------------------------------------------------------------------

void step_1()
{
//bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb
/*
	int i,j,sn;
	double Cy, Cy1, Cy2;

	//u - ��� ����� ��������� p1+p_k+ro*g*h; ks0 � ks1 - ��� k2*k/mu2 �������������� � ������ ����� � � ��������.
	struct 
	{
		double l, r, u, d;
	} u, ks0, ks1;

	//��������� �������������
	struct 
	{
		int l, r, u, d;
	} v;

						 //���������� ��������� ����:
	//���������� ���������� �����:
	for(i=1;i<N1-1;i++)
	{
		for(j=1;j<N2-1;j++)
		{
			//����� ����� � ������ ����� (�� �������� ����):
			sn=pnt[i][j].f_value[0].sand;

			//���������� v:
			v.l=pnt[i][j].vrt.l0;
			v.r=pnt[i][j].vrt.r0;
			v.u=pnt[i][j].vrt.u0;
			v.d=pnt[i][j].vrt.d0;

			//���������� u:
			u.r=-((p[i+1][j].value+pnt[i+1][j].f_value[v.r].p_k-p[i][j].value-pnt[i][j].f_value[v.r].p_k)/h1);
			u.l=-((p[i][j].value+pnt[i][j].f_value[v.l].p_k-p[i-1][j].value-pnt[i-1][j].f_value[v.l].p_k)/h1);
			u.d=-((p[i][j+1].value+pnt[i][j+1].f_value[v.d].p_k-p[i][j].value-pnt[i][j].f_value[v.d].p_k)/h2-ro2*g);
			u.u=-((p[i][j].value+pnt[i][j].f_value[v.u].p_k-p[i][j-1].value-pnt[i][j-1].f_value[v.u].p_k)/h2-ro2*g);

			//���������� ks0 � ks1
			ks0.l=pnt[i][j].f_value[v.l].k2*k[pnt[i][j].f_value[v.l].sand]/mu2;
			ks0.r=pnt[i][j].f_value[v.r].k2*k[pnt[i][j].f_value[v.r].sand]/mu2;
			ks0.u=pnt[i][j].f_value[v.u].k2*k[pnt[i][j].f_value[v.u].sand]/mu2;
			ks0.d=pnt[i][j].f_value[v.d].k2*k[pnt[i][j].f_value[v.d].sand]/mu2;
			ks1.l=pnt[i-1][j].f_value[v.l].k2*k[pnt[i-1][j].f_value[v.l].sand]/mu2;
			ks1.r=pnt[i+1][j].f_value[v.r].k2*k[pnt[i+1][j].f_value[v.r].sand]/mu2;
			ks1.u=pnt[i][j-1].f_value[v.u].k2*k[pnt[i][j-1].f_value[v.u].sand]/mu2;
			ks1.d=pnt[i][j+1].f_value[v.d].k2*k[pnt[i][j+1].f_value[v.d].sand]/mu2;

			//���������� Cy
			Cy1=((u.r+absd(u.r))*ks0.r-(u.l-absd(u.l))*ks0.l-(u.l+absd(u.l))*ks1.l+(u.r-absd(u.r))*ks1.r)/(2*h1);
			Cy2=((u.d+absd(u.d))*ks0.d-(u.u-absd(u.u))*ks0.u-(u.u+absd(u.u))*ks1.u+(u.d-absd(u.d))*ks1.d)/(2*h2);
			Cy=Cy1+Cy2;

			//���������� ���������� ���� �� �������:
			S_next[i][j].value[0]=S[i][j].value[0]-tau/m[sn]*Cy;

			//���� ���������, ��� �����:
			if(is_ist(i,j)) S_next[i][j].value[0]+=tau/m[sn]*Q/(h2*h2);


			//!!!!!!!!!!!!!!
			//������������� ������������� ��� ������������ ��������:
			if(Zaglushka)
			{
				if(S_next[i][j].value[0]>=1-S_wr[sn]-0.01) 
				{
					S_next[i][j].value[0]=1-S_wr[sn]-0.01;
					//printf("��������, ������������ ��������!\n");
					//getchar();
				}
			}
			//!!!!!!!!!!!!!!

		}
	}
*/
//eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee

	step_1_s_w();
	step_1_s_n();

	//���������� ��������� �����:
	step_1_gr();

	//���������� ������������ ����:
	//step_1_virt();

	step_1_s_g();

	//test_rec(2);
	//getchar();
}

//-------------------------------------------------------------------------------------------------------------------------------

void step_1_s_w()
{
	int i,j,sn;
	double Cy, Cy1, Cy2;

	//u - ��� ����� ��������� p1+p_k+ro*g*h; ks0 � ks1 - ��� k2*k/mu2 �������������� � ������ ����� � � ��������.
	struct 
	{
		double l, r, u, d;
	} u, ks0, ks1;

	//��������� �������������
	struct 
	{
		int l, r, u, d;
	} v;

	
	for(i=1;i<N1-1;i++)
	{
		for(j=1;j<N2-1;j++)
		{
			//����� ����� � ������ ����� (�� �������� ����):
			sn=pnt[i][j].f_value[0].sand;

			//���������� v:
			v.l=pnt[i][j].vrt.l0;
			v.r=pnt[i][j].vrt.r0;
			v.u=pnt[i][j].vrt.u0;
			v.d=pnt[i][j].vrt.d0;

			//���������� u:
			u.r=-((p[i+1][j].w[v.r]-p[i][j].w[v.r])/h1);
			u.l=-((p[i][j].w[v.l]-p[i-1][j].w[v.l])/h1);
			u.d=-((p[i][j+1].w[v.d]-p[i][j].w[v.d])/h2-ro_w*g);
			u.u=-((p[i][j].w[v.u]-p[i][j-1].w[v.u])/h2-ro_w*g);


			
			//���������� ks0 � ks1
			ks0.l=pnt[i][j].f_value[v.l].k_rw*k[pnt[i][j].f_value[v.l].sand]/mu_w;
			ks0.r=pnt[i][j].f_value[v.r].k_rw*k[pnt[i][j].f_value[v.r].sand]/mu_w;
			ks0.u=pnt[i][j].f_value[v.u].k_rw*k[pnt[i][j].f_value[v.u].sand]/mu_w;
			ks0.d=pnt[i][j].f_value[v.d].k_rw*k[pnt[i][j].f_value[v.d].sand]/mu_w;
			ks1.l=pnt[i-1][j].f_value[v.l].k_rw*k[pnt[i-1][j].f_value[v.l].sand]/mu_w;
			ks1.r=pnt[i+1][j].f_value[v.r].k_rw*k[pnt[i+1][j].f_value[v.r].sand]/mu_w;
			ks1.u=pnt[i][j-1].f_value[v.u].k_rw*k[pnt[i][j-1].f_value[v.u].sand]/mu_w;
			ks1.d=pnt[i][j+1].f_value[v.d].k_rw*k[pnt[i][j+1].f_value[v.d].sand]/mu_w;

			//���������� Cy
			Cy1=((u.r+absd(u.r))*ks0.r-(u.l-absd(u.l))*ks0.l-(u.l+absd(u.l))*ks1.l+(u.r-absd(u.r))*ks1.r)/(2*h1);
			Cy2=((u.d+absd(u.d))*ks0.d-(u.u-absd(u.u))*ks0.u-(u.u+absd(u.u))*ks1.u+(u.d-absd(u.d))*ks1.d)/(2*h2);
			Cy=Cy1+Cy2;

			//���������� ���������� ���� �� �������:
			S_next[i][j].w[0]=S[i][j].w[0]-tau/m[sn]*Cy;

			//�����:
			if(j==1 && Rain) S_next[i][j].w[0]+=tau/m[sn]*Q/(h2*h2);

/*
			//���� ���������, ��� �����:
			if(is_ist(i,j)) S_next[i][j].w[0]+=tau/m[sn]*Q/(h2*h2);


			//!!!!!!!!!!!!!!
			//������������� ������������� ��� ������������ ��������:
			if(Zaglushka)
			{
				if(S_next[i][j].value[0]>=1-S_wr[sn]-0.01) 
				{
					S_next[i][j].value[0]=1-S_wr[sn]-0.01;
					//printf("��������, ������������ ��������!\n");
					//getchar();
				}
			}
			//!!!!!!!!!!!!!!
*/
		}
	}

}

//-------------------------------------------------------------------------------------------------------------------------------

void step_1_s_n()
{
	int i,j,sn;
	double Cy, Cy1, Cy2;

	//u - ��� ����� ��������� p1+p_k+ro*g*h; ks0 � ks1 - ��� k2*k/mu2 �������������� � ������ ����� � � ��������.
	struct 
	{
		double l, r, u, d;
	} u, ks0, ks1;

	//��������� �������������
	struct 
	{
		int l, r, u, d;
	} v;

	
	for(i=1;i<N1-1;i++)
	{
		for(j=1;j<N2-1;j++)
		{
			//����� ����� � ������ ����� (�� �������� ����):
			sn=pnt[i][j].f_value[0].sand;

			//���������� v:
			v.l=pnt[i][j].vrt.l0;
			v.r=pnt[i][j].vrt.r0;
			v.u=pnt[i][j].vrt.u0;
			v.d=pnt[i][j].vrt.d0;

			//���������� u:
			u.r=-((p[i+1][j].n[v.r]-p[i][j].n[v.r])/h1);
			u.l=-((p[i][j].n[v.l]-p[i-1][j].n[v.l])/h1);
			u.d=-((p[i][j+1].n[v.d]-p[i][j].n[v.d])/h2-ro_n*g);
			u.u=-((p[i][j].n[v.u]-p[i][j-1].n[v.u])/h2-ro_n*g);


			
			//���������� ks0 � ks1
			ks0.l=pnt[i][j].f_value[v.l].k_rn*k[pnt[i][j].f_value[v.l].sand]/mu_n;
			ks0.r=pnt[i][j].f_value[v.r].k_rn*k[pnt[i][j].f_value[v.r].sand]/mu_n;
			ks0.u=pnt[i][j].f_value[v.u].k_rn*k[pnt[i][j].f_value[v.u].sand]/mu_n;
			ks0.d=pnt[i][j].f_value[v.d].k_rn*k[pnt[i][j].f_value[v.d].sand]/mu_n;
			ks1.l=pnt[i-1][j].f_value[v.l].k_rn*k[pnt[i-1][j].f_value[v.l].sand]/mu_n;
			ks1.r=pnt[i+1][j].f_value[v.r].k_rn*k[pnt[i+1][j].f_value[v.r].sand]/mu_n;
			ks1.u=pnt[i][j-1].f_value[v.u].k_rn*k[pnt[i][j-1].f_value[v.u].sand]/mu_n;
			ks1.d=pnt[i][j+1].f_value[v.d].k_rn*k[pnt[i][j+1].f_value[v.d].sand]/mu_n;

			//���������� Cy
			Cy1=((u.r+absd(u.r))*ks0.r-(u.l-absd(u.l))*ks0.l-(u.l+absd(u.l))*ks1.l+(u.r-absd(u.r))*ks1.r)/(2*h1);
			Cy2=((u.d+absd(u.d))*ks0.d-(u.u-absd(u.u))*ks0.u-(u.u+absd(u.u))*ks1.u+(u.d-absd(u.d))*ks1.d)/(2*h2);
			Cy=Cy1+Cy2;

			//���������� ���������� ���� �� �������:
			S_next[i][j].n[0]=S[i][j].n[0]-tau/m[sn]*Cy;


			//�������� ����������:
		/*	if(i==25 && j==8)
			{
				printf("u.d = %6.4e\n", u.d);
				printf("ks0.d = %6.4e\n", ks0.d);
				printf("ks1.u = %6.4e\n", ks1.u);
			}*/
/*
			//���� ���������, ��� �����:
			if(is_ist(i,j)) S_next[i][j].w[0]+=tau/m[sn]*Q/(h2*h2);


			//!!!!!!!!!!!!!!
			//������������� ������������� ��� ������������ ��������:
			if(Zaglushka)
			{
				if(S_next[i][j].value[0]>=1-S_wr[sn]-0.01) 
				{
					S_next[i][j].value[0]=1-S_wr[sn]-0.01;
					//printf("��������, ������������ ��������!\n");
					//getchar();
				}
			}
			//!!!!!!!!!!!!!!
*/
		}
	}

}

//-------------------------------------------------------------------------------------------------------------------------------

void step_1_s_g()
{
	int i,j;
	
	for(i=0;i<N1;i++)
	{
		for(j=0;j<N2;j++)
		{
			S_next[i][j].g[0] = 1 - S_next[i][j].w[0] - S_next[i][j].n[0];

			if(S_next[i][j].g[0]<0)
			{
				printf("S_g<0 v tochke i=%i j=%i\n",i,j);
				getchar();
			}

		}
	}

}

//-------------------------------------------------------------------------------------------------------------------------------


void step_1_gr()
{
	int i, j;

	//������� �������:
	for(i=0;i<N1;i++)	
	{	
		S_next[i][0].w[0]=S_next[i][1].w[0];
		S_next[i][0].n[0]=S_next[i][1].n[0];		
		if(Gr_usl && is_ist(i,1)) S_next[i][0].n[0]=S_n_ist;
	}

	//������� �������
	for(j=0;j<N2;j++)	
	{
		S_next[0][j].n[0]=S_next[1][j].n[0];
		S_next[N1-1][j].n[0]=S_next[N1-2][j].n[0];

		S_next[0][j].w[0]=S_next[1][j].w[0];
		S_next[N1-1][j].w[0]=S_next[N1-2][j].w[0];
	}
	//������ �������
	j=N2-1;
	for(i=0;i<N1;i++)	
	{
		
		S_next[i][j].n[0]=S_next[i][j-1].n[0];
		S_next[i][j].w[0]=S_next[i][j-1].w[0];
	}	


	//step_1_gr_pres();
	
}

//-------------------------------------------------------------------------------------------------------------------------------

void step_1_gr_pres()
{
	int i, j;

	//������� �������:
	for(i=0;i<N1;i++)	
	{	
		p[i][0].w[0] = p[i][1].w[0] - ro_w*g;
		if(!(Gr_usl && is_ist(i,1))) p[i][0].n[0] = p[i][1].n[0] - ro_n*g;		
		//if(Gr_usl && is_ist(i,1)) S_next[i][0].n[0]=S_n_ist;
	}

	//������� �������
	for(j=0;j<N2;j++)	
	{
		p[0][j].n[0]=p[1][j].n[0];
		p[N1-1][j].n[0]=p[N1-2][j].n[0];

		p[0][j].w[0]=p[1][j].w[0];
		p[N1-1][j].w[0]=p[N1-2][j].w[0];
	}
	//������ �������
	/*j=N2-1;
	for(i=0;i<N1;i++)	
	{
		
		p[i][j].n[0]=p[i][j-1].n[0] + ro_n*g;
		p[i][j].w[0]=p[i][j-1].w[0] + ro_w*g;
	}*/	
	
}

//-------------------------------------------------------------------------------------------------------------------------------

void step_1_virt()
{
	//�-�� ������������ ������ ����� ���������� ��������� ����

	int i, j;
	double temp;

	for(i=0;i<N1;i++)
	{
		for(j=0;j<N2;j++)
		{
			//��� ���� ����� ������� �������� � ��������� ����:
			S_next[i][j].value[1]=S_next[i][j].value[0];

			//��� ����� �� ���������� ����������� ������������ �������:
			if(is_interf(i,j) && Interf_usl)
			{
				if(S_next[i][j].value[0]<1-S_zv)	S_next[i][j].value[1]=0;
				else
				{
					//
					//printf("interf %i %i \n",i,j);
					if(1-S_next[i][j].value[0]-S_wr[0]<=0)
					{
						printf("\n\n\n\nAchtung!!!\n\n\n");
						getchar();
					}
					temp=powf(p_d[1]/p_d[0],lambda[1])*powf((1-S_next[i][j].value[0]-S_wr[0])/(1-S_wr[0]),lambda[1]/lambda[0]);
					S_next[i][j].value[1]=(1-S_wr[1])*(1-temp);
				}
			}
		}
	}

}

//-------------------------------------------------------------------------------------------------------------------------------

int is_interf(int i,int j)
{
	if(pnt[i][j].f_value[1].sand!=pnt[i][j].f_value[0].sand) return 1;
	else return 0;
}

//-------------------------------------------------------------------------------------------------------------------------------

int is_ist(int i, int j)
{
	if(j==1 && abs(i-N1/2)<=N_ist/2) return 1;
	else return 0;
}

//-------------------------------------------------------------------------------------------------------------------------------

/*void step_2()
{
	int ii;
	//��������� � ������� �������� �������:
	double nev0,nev;
	//���� ���������� ��������� ������� �� ������ ��������:
	const int flag_pr_nev=0;

	//��������� ������� ������� ����� ���������� ������ ��������, � ���� �������� 1
	nev0=1;
	nev=nev0;

	for(ii=0;nev/nev0>epsilon;ii++)
	{
		step_2_relax_iter();
		nev=step_2_nevazka();

		if(flag_pr_nev) printf("nev/nev0 = %6.4e\n",nev/nev0);

		//��������� ������� ������� ����� ���������� ������ ��������:
		if(ii==0)
		{
			nev0=nev;
			//���� ����� �������� ������� �������, �� ����� ����� �� �����:
			if(nev0==0) nev0=1;
			//���� ����� ���, �� �� ������ ����� �� ������� ������������ �������� �������:
			//if(Linz_exist==0 && vrema<2) nev0=1e-5;
		}
	}

	//����� �� ����� ������ ���� �� �������:
	printf("vrema = %i;\t",vrema);
	//����� ���-�� ��������
	printf("iter1 = %i;\n",ii);


}

//-------------------------------------------------------------------------------------------------------------------------------

void step_2_relax_iter()
{
	double chisl;
	int i,j;
	//�������� ��� ����������:
	const double omega=1.3;	

		
	//��������� ��������� ��������� �����:
	step_2_gr_before();

	//���������� �����
	for(i=1;i<N1-1;i++)	
	{
		for(j=1;j<N2-1;j++)
		{
			
			chisl=(pnt[i][j].rel.pravo*p[i+1][j].value+pnt[i][j].rel.levo*p_next[i-1][j].value+pnt[i][j].rel.niz*p[i][j+1].value+pnt[i][j].rel.verh*p_next[i][j-1].value+pnt[i][j].rel.f)*omega+pnt[i][j].rel.znam*p[i][j].value*(1-omega);
			p_next[i][j].value=chisl/pnt[i][j].rel.znam;
		}
	}

	//������ ��������� �������:
	step_2_gr_after();
	//step_2_gr_after_2();


	//������ ������� p � p_next:
	p_buf=p;
	p=p_next;
	p_next=p_buf;
		
}

//-------------------------------------------------------------------------------------------------------------------------------

void step_2_gr_before()
{
	int i, j;

	//������� ���� ��������� � ����������� ����
	for(i=0;i<N1;i++)
		p_next[i][0].value=p[i][0].value;

	//����� ���� ��������� � ����������� ����
	for(j=0;j<N2;j++)	
		p_next[0][j].value=p[0][j].value;

}

//-------------------------------------------------------------------------------------------------------------------------------

void step_2_gr_after()
{
	int i, j;

	for(j=0;j<N2;j++)
	{
		p_next[0][j].value=p_atm+ro1*g*j*h2;
		p_next[N1-1][j].value=p_atm+ro1*g*j*h2;
	}

	j=N2-1;	//������ � ������� �������
	for(i=0;i<N1;i++)
	{
		//if(is_ist(i,1)) p_next[i][0].value=p_next[i][1].value;//-h2*ro1*g;
		//else p_next[i][0].value=p_atm;
		p_next[i][0].value=p_next[i][1].value;//-h2*ro1*g;
		p_next[i][j].value=p_next[i][j-1].value;//+h2*ro1*g;
	}
}

//-------------------------------------------------------------------------------------------------------------------------------

void step_2_gr_after_2()
{
	int i, j;

	i=N1-1;
	for(j=0;j<N2;j++)
	{
		p_next[0][j].value=p_next[1][j].value-(p_next[2][j].value-p_next[1][j].value)*(pnt[1][j].f_value[0].k1+pnt[2][j].f_value[0].k1)/(pnt[0][j].f_value[0].k1+pnt[1][j].f_value[0].k1);
		p_next[i][j].value=p_next[i-1][j].value-(p_next[i-2][j].value-p_next[i-1][j].value)*(pnt[i-1][j].f_value[0].k1+pnt[i-2][j].f_value[0].k1)/(pnt[i][j].f_value[0].k1+pnt[i-1][j].f_value[0].k1);
	}

	j=N2-1;	//������ � ������� �������
	for(i=0;i<N1;i++)
	{
		//if(is_ist(i,1)) p_next[i][0].value=p_next[i][1].value;//-h2*ro1*g;
		//else p_next[i][0].value=p_atm;
		p_next[i][0].value=p_next[i][1].value-h2*ro1*g;
		p_next[i][j].value=p_next[i][j-1].value+h2*ro1*g;
	}
}


//-------------------------------------------------------------------------------------------------------------------------------

double step_2_nevazka()
{
	double chisl;
	int i,j;
	double value=0;

	//������� ������� ������ �� ���������� ������:
	//���������� �����
	for(i=1;i<N1-1;i++)	
	{
		for(j=1;j<N2-1;j++)
		{
			
			chisl=(pnt[i][j].rel.pravo*p[i+1][j].value+pnt[i][j].rel.levo*p[i-1][j].value+pnt[i][j].rel.niz*p[i][j+1].value+pnt[i][j].rel.verh*p[i][j-1].value+pnt[i][j].rel.f);
			value+=(pnt[i][j].rel.znam*p[i][j].value-chisl)*(pnt[i][j].rel.znam*p[i][j].value-chisl)*h1*h2;
			
		}
	}

	value=sqrt(value);
	return value;
}*/

//-------------------------------------------------------------------------------------------------------------------------------

void save()
{
	FILE *fp;
	int numb,num,i;
	char namefl[15];	
	const int vrs=Version;

	numb=vrema/M_out;
	sprintf(namefl, "save\\S%d.set",numb+1000);
	if((fp=fopen(namefl,"wb"))==NULL)
		printf("�� ������\n\n\n\n");

	num=0;
	//������ ������ ������:
	num += fwrite( &vrs, sizeof( int ), 1, fp );
	//������ ������������:
	for(i=0;i<N1;i++) num +=  fwrite( S[i], sizeof( Sat ), (N2), fp );
	//������ ��������:
	for(i=0;i<N1;i++) num +=  fwrite( p[i], sizeof( Pres ), (N2), fp );
	//������ ������ ���� �� �������:
	num += fwrite( &vrema, sizeof( int ), 1, fp );

	if(num!=N1*(N2)+N1*(N2)+2) printf("Oshibka sohranenia\n\n\n\n\n\n\n\n");

	fclose(fp);

}

//-------------------------------------------------------------------------------------------------------------------------------

void load()
{
	FILE *fp;
	int num,i;
	int vrs;

	if((fp=fopen("save\\S1119.set","rb"))==NULL)
		printf("�� ������\n\n\n\n");

	num=0;

	num += fread( &vrs, sizeof( int ), 1, fp );
	for(i=0;i<N1;i++) num +=  fread( S[i], sizeof( Sat ), (N2), fp );
	for(i=0;i<N1;i++) num +=  fread( p[i], sizeof( Pres ), (N2), fp );
	num += fread( &vrema, sizeof( int ), 1, fp );

	if(num!=N1*(N2)+N1*(N2)+2) printf("Oshibka zagruzki\n\n\n\n\n\n\n\n");

	fclose(fp);

	vrema++;
	
}

//-------------------------------------------------------------------------------------------------------------------------------

void rec_info()
{
	
	//��������� �� ����:
	FILE *f1;

	
	//������ ��� ����� �����:
	char namefl[15];
	
	
					//��������� ��� ����� � ������� ���, � ��������� � ����������:
	sprintf(namefl, "out_test\\info.dat");
    f1=fopen(namefl,"w");

					//������ � ����:	
	
	fprintf(f1,"�������� ����� �� ������������:\n");
	fprintf(f1,"h1 = %6.4e\n",h1);
	fprintf(f1,"h2 = %6.4e\n",h2);

	fprintf(f1,"���������� (������ ����� � �����) ������� ����� �����, �.�. ������� � ������:\n");
	fprintf(f1,"N_x_1 = %i\n",N_x_1);
	fprintf(f1,"N_x_2 = %i\n",N_x_2);
	fprintf(f1,"N_y_1 = %i\n",N_y_1);
	fprintf(f1,"N_y_2 = %i\n",N_y_2);
	

	//�������� ������:
	fclose(f1);
	
}










//-------------------------------------------------------------------------------------------------------------------------------


void test_print(int pr_type)
{
	//�������� pr_type ����������, ��� ����� ���������������. ���� �� ������ 100, �� ������. �������� �� ���� ������.
	int i,j;


	if(pr_type<100)
	{
			//���������� ����������� ��������� S �� �������:
		//printf("%f\n\n\n",1-S_zv);

		//������ ������� N_x_1, N_x_2, N_y_1, N_y_2
		printf("%f\n",x[N_x_1]);
		printf("%f\n",x[N_x_2]);
		printf("%f\n",y[N_y_1]);
		printf("%f\n",y[N_y_2]);
	}

	if(pr_type>=100)
	{
		for(j=0;j<N2;j++)
		{
			for(i=0;i<N1;i++)
			{
				//���� ����� �� ���� ������ � �� �����������
				if(pr_type==101) printf("%i%i ", pnt[i][j].cl_x, pnt[i][j].cl_y);
				if(pr_type==102) printf("%i ", pnt[i][j].cl_x);
				if(pr_type==103) printf("%i ", pnt[i][j].cl_y);
			}
			printf("\n");
		}
	}

}

//-------------------------------------------------------------------------------------------------------------------------------

void test_rec(int rec_type)
{
	//�������� rec_type ����������, ��� ����� ���������������. 
	//���� ������ 100, �� ��� ��������� ����������� 
	//��� ����� ����� ������ �������� �� ��������� vrema

	//��������� �� ����:
	FILE *f1;

	int numb,j,i;
	//������ ��� ����� �����:
	char namefl[15];
	//���������� ��� �������� ��������:
	double value1=0;
	
	numb=vrema/M_out;

					//��������� ����� ������ � ������� ��, � ��������� � �����������:
	sprintf(namefl, "out_test\\t_%d_%d.dat",rec_type,vrema);
    f1=fopen(namefl,"w");

					//������ � �����:	
	//������ �����:
	fprintf(f1,"TITLE =  \"Test\" \n");
	fprintf(f1,"VARIABLES = \"X\",\"Y\",\"Test\" \n");
    fprintf(f1,"ZONE T = \"BIG ZONE\", J=%d  ,I=%d, F = POINT \n", N2, N1);
	
	//������ ������:
	for(j=0; j<N2; j++)
	{
			for(i=0; i<N1; i++)
			{
				if(rec_type==1) value1=S[i][j].value[0];
				if(rec_type==2) value1=S[i][j].value[1];
				if(rec_type==3) value1=S_next[i][j].value[0];
				if(rec_type==4) value1=S_next[i][j].value[1];
				if(rec_type==5) value1=p[i][j].value;
				if(rec_type==6) value1=p_next[i][j].value;
				if(rec_type==7) value1=pnt[i][j].cl_x;
				if(rec_type==8) value1=pnt[i][j].cl_y;
				if(rec_type==9) value1=pnt[i][j].vrt.l0;
				if(rec_type==10) value1=pnt[i][j].vrt.r0;
				if(rec_type==11) value1=pnt[i][j].vrt.u0;
				if(rec_type==12) value1=pnt[i][j].vrt.d0;
				if(rec_type==13) value1=pnt[i][j].f_value[0].sand;
				if(rec_type==14) value1=pnt[i][j].f_value[1].sand;
				if(rec_type==15) value1=is_interf(i, j);

			     fprintf(f1,"%6.4e %6.4e %6.4e \n",x[i],y[j],value1);
				 
			}
	 }


	//�������� ���� �������� ������:
	fclose(f1);
	
}

//-------------------------------------------------------------------------------------------------------------------------------

void test1()
{
	int i;
	for(i=1;i<=15;i++)
		test_rec(i);
}

//-------------------------------------------------------------------------------------------------------------------------------

double absd(double val)
{
	if(val<0) return -val;
	else return val;
}

//-------------------------------------------------------------------------------------------------------------------------------

void test_step_1()
{
	int i,j,sn;
	double Cy, Cy1, Cy2;

	//u - ��� ����� ��������� p1+p_k+ro*g*h; ks0 � ks1 - ��� k2*k/mu2 �������������� � ������ ����� � � ��������.
	struct 
	{
		double l, r, u, d;
	} u, ks0, ks1;

	//��������� �������������
	struct 
	{
		int l, r, u, d;
	} v;

						 //���������� ��������� ����:
	//���������� ���������� �����:
	for(i=1;i<N1-1;i++)
	{
		for(j=1;j<N2-1;j++)
		{
			//����� ����� � ������ ����� (�� �������� ����):
			sn=pnt[i][j].f_value[0].sand;

			//���������� v:
			v.l=pnt[i][j].vrt.l0;
			v.r=pnt[i][j].vrt.r0;
			v.u=pnt[i][j].vrt.u0;
			v.d=pnt[i][j].vrt.d0;

			//���������� u:
			u.r=-((p[i+1][j].value+pnt[i+1][j].f_value[v.r].p_k-p[i][j].value-pnt[i][j].f_value[v.r].p_k)/h2);
			u.l=-((p[i][j].value+pnt[i][j].f_value[v.l].p_k-p[i-1][j].value-pnt[i-1][j].f_value[v.l].p_k)/h2);
			u.d=-((p[i][j+1].value+pnt[i][j+1].f_value[v.d].p_k-p[i][j].value-pnt[i][j].f_value[v.d].p_k)/h2-ro2*g);
			u.u=-((p[i][j].value+pnt[i][j].f_value[v.u].p_k-p[i][j-1].value-pnt[i][j-1].f_value[v.u].p_k)/h2-ro2*g);

			//���������� ks0 � ks1
			ks0.l=pnt[i][j].f_value[v.l].k2*k[pnt[i][j].f_value[v.l].sand]/mu2;
			ks0.r=pnt[i][j].f_value[v.r].k2*k[pnt[i][j].f_value[v.r].sand]/mu2;
			ks0.u=pnt[i][j].f_value[v.u].k2*k[pnt[i][j].f_value[v.u].sand]/mu2;
			ks0.d=pnt[i][j].f_value[v.d].k2*k[pnt[i][j].f_value[v.d].sand]/mu2;
			ks1.l=pnt[i-1][j].f_value[v.l].k2*k[pnt[i-1][j].f_value[v.l].sand]/mu2;
			ks1.r=pnt[i+1][j].f_value[v.r].k2*k[pnt[i+1][j].f_value[v.r].sand]/mu2;
			ks1.u=pnt[i][j-1].f_value[v.u].k2*k[pnt[i][j-1].f_value[v.u].sand]/mu2;
			ks1.d=pnt[i][j+1].f_value[v.d].k2*k[pnt[i][j+1].f_value[v.d].sand]/mu2;

			//���������� Cy
			Cy1=((u.r)*ks0.r-(u.l)*ks0.l-(u.l)*ks1.l+(u.r)*ks1.r)/(2*h1);
			Cy2=((u.d)*ks0.d-(u.u)*ks0.u-(u.u)*ks1.u+(u.d)*ks1.d)/(2*h2);
			Cy=Cy1+Cy2;

			//���������� ���������� ���� �� �������:
			S_next[i][j].value[0]=S[i][j].value[0]-tau/m[sn]*Cy;

			//���� ���������, ��� �����:
			if(is_ist(i,j)) S_next[i][j].value[0]+=tau/m[sn]*Q/(h2*h2);


			//!!!!!!!!!!!!!!
			//������������� ������������� ��� ������������ ��������:
			if(S_next[i][j].value[0]>=1-S_wr[1]) S_next[i][j].value[0]=1-S_wr[1]-0.001;
			//!!!!!!!!!!!!!!

		}
	}

	//���������� ��������� �����:
	step_1_gr();

	//���������� ������������ ����:
	step_1_virt();
}






//-------------------------------------------------------------------------------------------------------------------------------

void bu_rec()
{
	//��������� �� �����:
	FILE *fS2,*fp1,*fp2,*f3S, *f3p;

	int numb,j,i;
	//������ ��� ����� �����:
	char namefl[15];	
	numb=vrema/M_out;

					//��������� ����� ������ � ������� ��, � ��������� � �����������:
	sprintf(namefl, "out_s\\S%d.dat",numb+1000);
	if((fS2=fopen(namefl,"w"))==NULL)
		printf("�� ������\n\n\n\n");
	sprintf(namefl, "out_p\\P%d.dat",numb+1000);
	if((fp1=fopen(namefl,"w"))==NULL)
		printf("�� ������\n\n\n\n");
	//��� �������� �����
	sprintf(namefl, "out_pn\\P%d.dat",numb+1000);
	if((fp2=fopen(namefl,"w"))==NULL)
		printf("�� ������\n\n\n\n");

	//��� 3d �����������:
	sprintf(namefl, "out1s\\S2_%d.dat",vrema);
    f3S=fopen(namefl,"w");
	sprintf(namefl, "out1p\\p1_%d.dat",vrema);
    f3p=fopen(namefl,"w");

					//������ � �����:
	//������� ���������:
	for(j=0;j<N2;j++)
	{
		//l=sand(j);
		fprintf(fS2,"%6.4e ",y[j]);
		fprintf(fS2,"%8.6e\n",S[N2/2][j].value[0]);
		fprintf(fp1,"%6.4e ",y[j]);
		fprintf(fp1,"%8.6e\n",p[N2/2][j].value);
		//fprintf(fp2,"%6.4e ",y[j]);
		//fprintf(fp2,"%8.6e\n",p1_a[1][j]+re[1][j].p_k);
	}

	//������ ����������:
	//������ �����:
	fprintf(f3S,"TITLE =  \"Saturation\" \n");
	fprintf(f3S,"VARIABLES = \"X\",\"Y\",\"S2\" \n");
    fprintf(f3S,"ZONE T = \"BIG ZONE\", J=%d  ,I=%d, F = POINT \n", N2, N1);

	fprintf(f3p,"TITLE =  \"Pressure\" \n");
	fprintf(f3p,"VARIABLES = \"X\",\"Y\",\"P1\" \n");
    fprintf(f3p,"ZONE T = \"BIG ZONE\", J=%d  ,I=%d, F = POINT \n", N2, N1);

	//������ ������:
	for(j=0; j<N2; j++)
	{
			for(i=0; i<N1; i++)
			{
			     fprintf(f3S,"%6.4e %6.4e %6.4e \n",x[i],y[j],S[i][j].value[0]);
				 fprintf(f3p,"%6.4e %6.4e %6.4e \n",x[i],y[j],p[i][j].value);
			}
	 }


	//�������� ���� �������� ������:
	fclose(fS2);
	fclose(fp1);
	fclose(fp2);
	fclose(f3S);
	fclose(f3p);

}

