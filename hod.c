#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_sf_expint.h>
#include <unistd.h>
#include "mpi.h"

//#define k_num  237
//#define k_num  10001
#define r_num  451
#define cut_number 6
//#define velocity_num 442
#define node   16
#define pi     3.141592653589793
#define G   6.67259E-11
#define little_h     0.7
#define H0    100*little_h
#define delta_c   1.686
//#define z   1.75
//#define Omega_m0   0.3
//#define Omega_lambda0   0.7
#define triangle   200.0
#define sqrt_2pi   2.506628275
#define parallel_max   500
#define perpendicular_max   28
//#define perpendicular_min   15
//#define number_density  18807/(80284132361.29527*16.4992*87.0/100.0/41252.96125)
//#define number_density  0.0006732289692270188
//#define sigma_ng  0.00006732289692270188
#define number  3
#define node_pho 100
#define node_vel 1000
#define sigma_log_M 0.2
#define alpha 1.0
#define m_max   70          //compare to m_max=70000,ka_square 30.0476 and 30.0458 ,but when the largest m is 15,ka_square is about 38,big difference
#define ndim  3
//#define walks  50000
#define walks 1
#define chains 1
#define dpi_max 100
//#define r_test 431



char Linear_matter_power_spectrum[200];     
char Nonlinear_matter_power_spectrum[200]; 
int  Start_point;
int  Perpendicular_min;	
char Wp_mean[200];
char Inverse_matrix[200];
double Number_density;			
double Sigma_ng;	
double Redshift;
int K_num;
double Omega_m0;
double Omega_lambda0;

double M[m_max]={0},R[m_max]={0},M_add[m_max+1]={0},R_add[m_max+1]={0};
double Sigma[m_max]={0},Sigma_add[m_max+1]={0},F[m_max]={0},Dlninvsigma_divide_dm[m_max]={0},Dn_divide_dm[m_max]={0},Bias[m_max]={0};
double Xi[node]={0},Wi[node]={0};
double R_corr[r_num]={0};
double Velxi[node_vel]={0},Velwi[node_vel]={0};
//double Velocity[velocity_num]={0};
double *K, *Pkl, *Pk;
double Parallel[parallel_max]={0}, *Perpendicular;
double Rh=0;
int r_cut[cut_number+1]= {0, 50, 80, 150, 200, 300, r_num};
//int r_cut_value[cut_number]={15, 20, 25, 30, 40, 50, 60, 80, 100, 100, 500, 1000, 5000, 10000, 100000};
int r_cut_value[cut_number]={100, 500, 1000, 3000, 3000, 5000};

struct Linear
{		
	double kl_r;
	double pkl_r;
};
		
struct Nonlinear
{		
	double k_r;
	double pk_r;
};



double Rho_m(double redshift)
{

	double E,Omega_m,Omega_lambda,H;
	double rho_m;
	E = sqrt(Omega_lambda0+Omega_m0*pow(1+redshift,3));
//	WRITE(*,*)Omega_m0,Omega_lambda0
//	printf(*,*)'E=',E
    Omega_m = Omega_m0*pow((1+redshift),3)/pow(E,2);
//	 printf(*,*)'Omega_m=',Omega_m
    Omega_lambda = Omega_lambda0/pow(E,2);
//	printf(*,*)'Omega_lambda=',Omega_lambda
    H = H0*E*1000/(3.08568E22);
//	printf(*,*)'H=',H
    rho_m = 3.0*H*H/(8.0*pi*G);
 // printf(*,*)'Rho_m=',Rho_m
    return  rho_m;
}



void Sigma_f(double *sigma_s, double *R, int length)
{
    int i,j,num;
	double xi_res[node],wi_res[node];
	gsl_interp_accel *acc0= gsl_interp_accel_alloc ();
	gsl_spline *spline0= gsl_spline_alloc (gsl_interp_cspline, K_num);
	gsl_spline_init (spline0, K, Pkl, K_num);
	gsl_integration_glfixed_table *table = gsl_integration_glfixed_table_alloc(node);

	for(num=0; num<node; num++)
	{
		gsl_integration_glfixed_point(0,2*pi,num,&Xi[num],&Wi[num],table);
//		printf("%lf\t%lf\n", xi[num], wi[num]);
	}
	for(j=0; j<length; j++)
	{
		sigma_s[j] = 0;
		for(i=0; i<100; i++)
		{	
			for(num=0; num<node; num++)
			{
				xi_res[num] = Xi[num] + 2*pi*i;
				wi_res[num] = Wi[num];
				sigma_s[j] = sigma_s[j] + wi_res[num]*xi_res[num]*xi_res[num]/R[j]/R[j]/R[j]*\
				pow(3.0*(sin(xi_res[num])-xi_res[num]*cos(xi_res[num]))/pow(xi_res[num],3),2)*gsl_spline_eval(spline0, xi_res[num]/R[j], acc0);
			}
		} 
		sigma_s[j] = sqrt(sigma_s[j]/2.0/pi/pi);
	}

}

void Multiplicity(double *f)
{
	int i;
	double A,a_temp,beta,b,c;
//	A = (0.1*log10(triangle)-0.05)*pow(1+Redshift,-0.14);
	A = 0.186*pow(1+Redshift,-0.14);
//	printf(*,*)"A=",A;
//	a_temp = (1.43+pow(log10(triangle)-2.3,1.5))*pow(1+Redshift,-0.06);
	a_temp = 1.47*pow(1+Redshift,-0.06);
// 	printf(*,*)"a_temp=",a_temp;
    beta = pow(10,-pow(0.75/log10(triangle/75.0),1.2));
//	printf(*,*)"beta=",beta;
//	b = (1.0+pow(log10(triangle)-1.6,-1.5))*pow(1+Redshift,-beta);
	b = 2.57*pow(1+Redshift,-beta);
//	printf(*,*)"b=",b;
//	c = 1.2+pow(log10(triangle)-2.35,1.6);
	c = 1.19;
//	printf(*,*)"c=",c;
    for(i=0; i<m_max; i++)
    {
    	f[i] = A*(pow(Sigma[i]/b,-a_temp)+1)*exp(-c/pow(Sigma[i],2));
    }
}

void bias_f(double *bias)
{
	int i;
	double nu[m_max];
	double y, A, a_temp, B, b_temp, C, c_temp;
    y = log10(triangle);
    A = 1 + 0.24*y*exp(-pow(4.0/y,4));
    a_temp = 0.44*y - 0.88;
    B = 0.183;
    b_temp = 1.5;
    C = 0.019 + 0.107*y +0.19*exp(-pow(4.0/y,4));
    c_temp = 2.4;
	for(i=0; i<m_max; i++)
	{
		nu[i] = delta_c/Sigma[i];      
		bias[i] = 1 - A*pow(nu[i],a_temp)/(pow(nu[i],a_temp)+pow(delta_c,a_temp)) + B*pow(nu[i],b_temp) + C*pow(nu[i],c_temp);
    }
}

void profile(double *prof, double *rvir)
{
	int i,j,s;
	double m_nl = 3.9E12;
    double c[m_max], rho_h, rs[m_max], rhos[m_max];
    double *si_more, *ci_more, *si, *ci;
   	si_more = (double *)malloc(m_max*K_num*sizeof(double));
   	ci_more = (double *)malloc(m_max*K_num*sizeof(double));
   	si = (double *)malloc(m_max*K_num*sizeof(double));
   	ci = (double *)malloc(m_max*K_num*sizeof(double));
	rho_h = triangle*Rh;
    for(s=0; s<m_max; s++)
    {
		rvir[s] = pow(3.0*M[s]/4.0/pi/rho_h,1.0/3.0);
    }

	for(s=0; s<m_max; s++)
	{
//		printf("red%lf\n",Redshift[i]);
//		printf("m_nl%le\n", m_nl[i]);
//		printf("m%le\n", M[s]);
		c[s] = 11.0*pow(M[s]/m_nl,-0.13)/(1+Redshift);
		rs[s] = rvir[s]/c[s];
		rhos[s] = Rh*triangle/3.0*pow(c[s],3)/(log(1+c[s])-c[s]/(1+c[s]));
//		printf("rho_h = %lf\trvir= %lf\tc = %lf\trs = %lf\trhos = %lf\n ", rho_h, rvir[s],  c[i*m_max+s], rs[i*m_max+s], rhos[i*m_max+s]);
    	for(j=0; j<K_num; j++)
    	{
			si_more[s*K_num+j] = gsl_sf_Si((1+c[s])*K[j]*rs[s]);
			ci_more[s*K_num+j] = gsl_sf_Ci((1+c[s])*K[j]*rs[s]);               
//			printf(*,*)'si_more=',si_more
//			printf(*,*)'ci_more=',ci_more
			si[s*K_num+j] = gsl_sf_Si(K[j]*rs[s]);
			ci[s*K_num+j] = gsl_sf_Ci(K[j]*rs[s]);
//			printf(*,*)'si=',si
//			printf(*,*)'ci=',ci
			prof[s*K_num+j] = 4.0*pi*rhos[s]*pow(rs[s],3)/M[s]*( sin(K[j]*rs[s]) * (si_more[s*K_num+j]-si[s*K_num+j]) - 
			sin(c[s]*K[j]*rs[s])/(1+c[s])/K[j]/rs[s] + cos(K[j]*rs[s]) * (ci_more[s*K_num+j]-ci[s*K_num+j]) ); 
//			printf("prof[%d] = %.10lf\n", j, (*prof)[j]);
		}
	}
	free(si_more);
	free(ci_more);
	free(si);
	free(ci);
}

void exclusion(double *exc, double *rvir)
{
	int i,j;
	int loc;
	for(i=0; i<m_max; i++)
	{
		loc = (int)100.0*(log10(2.0*rvir[i])+1.005);
		for(j=0; j<r_num; j++)
		{
			if(j<=loc)
			{
				exc[i*r_num+j] = 0;
			}
			else 
			{
				exc[i*r_num+j] = 1;
			}
		}
	}
}

double approximation(double pa_in, double pe_in, double v_in)
{
		double distance;
		int loc;
		distance = sqrt(pow(pe_in,2)+pow(pa_in-v_in,2));
		loc = (int)(log10(distance)/0.01+100.5);
		if (loc>=r_num) loc = r_num-1;
//		if (loc>=r_test) loc = r_test-1;
		return loc;		
}







double getwp(double *param, double *w_analysis, double *halo_mass ,double *fsat, double *beff)
{
	int i, j, mum, num, pa, pe, loop;
	int location;
	double ng=0,nc=0,ns=0;
	double temp1[m_max]={0};
	double N_cen_average[m_max]={0},N_sat_average[m_max]={0},N_average[m_max]={0};
	double rvir[m_max]={0};
	double *u;
	double *p1;
	double *pcs1, *pss1, *p;
	double exc[m_max*r_num]={0},fexc[r_num]={0};
	double corr_mm[r_num]={0},corr_1h[r_num]={0},corr_2h[r_num]={0},corr[r_num]={0};
	double corr_2hold1[r_num]={0},corr_2hold2[r_num]={0},corr_2hnew[r_num]={0};
	double *xi_s[node_pho*number];
	double phoxi[node_pho*number]={0};	
	double phowi[node_pho*number]={0};
	double log_M_min, M0, M1_primer, error;
	double double_sigma_square;
	double xi_res[node]={0},wi_res[node]={0};
	double corrrrr[r_num]={0};
	double *w[number];
	double low_limit,high_limit;
	gsl_integration_glfixed_table *table_pho;

	p1 = (double *)malloc(K_num*sizeof(double));
	pcs1 = (double *)malloc(K_num*sizeof(double));
	pss1 = (double *)malloc(K_num*sizeof(double));
	p = (double *)malloc(K_num*sizeof(double));
	u = (double *)malloc(m_max*K_num*sizeof(double));
	for(i=0 ;i<node_pho*number; i++)
	{
		xi_s[i] = (double *)malloc(Perpendicular_min*sizeof(double));
	}
	for(i=0 ;i<node_pho*number; i++)
	{
		for(j=0; j<Perpendicular_min; j++)
		{
			xi_s[i][j] = 0;
		}
	}	
	for(i=0 ;i<number; i++)
	{
		w[i] = (double *)malloc(Perpendicular_min*sizeof(double));
	}	

	for(i=0 ;i<number; i++)
	{
		for(j=0; j<Perpendicular_min; j++)
		{
			w[i][j] = 0;
		}
	}

	for(i=0 ;i<K_num; i++)
	{
		p1[i] = 0;
		pcs1[i] = 0;
		pss1[i] = 0;
		p[i] = 0;
	}

//	FILE *ft;
//	double corrrrr[r_test]={0};

	log_M_min = *(param+0);
	M1_primer = pow( 10.0,*(param+1) );
	M0 		  = pow( 10.0,0.76*(*(param+1))+2.3 );
	error 	  =	*(param+2);
	double_sigma_square = error*error*2.0;
//	printf("log_M_min=%lf\tsigma_log_M=%lf\tM0=%le\tM1_primer=%le\talpha=%lf\n", log_M_min, sigma_log_M, M0, M1_primer, alpha);
//	printf("error=%lf\n",error);
	profile(u, rvir);
	exclusion(exc, rvir);	
	for(i=0;i<m_max;i++)
	{
		temp1[i] = (log10(M[i]) - log_M_min)/sigma_log_M;
//		printf(*,*)'temp1=',temp1
		N_cen_average[i] = (gsl_sf_erf(temp1[i])+1.0)/2.0;
//		printf(*,*)'N_cen_average=',N_cen_average
		if((M[i] - M0)<=0)
		{
			N_sat_average[i]=0;
		}
		else        
		{
			N_sat_average[i] = N_cen_average[i]*pow((M[i] - M0)/M1_primer,alpha); 		
		}       
//		printf(*,*)'N_sat_average=',N_sat_average
		N_average[i] = N_cen_average[i] + N_sat_average[i];
//		printf(*,*)'N_average=',N_average
//		for(j=0; j<k_num; j++)
//		{
//			printf("u[%d][%d]=%lf\n", i, j ,u[i][j]);
//		}
//		printf("N_cen_average[%d] = %le\tN_sat_average[%d] = %le\tN_average[%d] = %le\n", i, N_cen_average[i], i, N_sat_average[i], i, N_average[i]);
//		printf("rvir[%d] = %lf\n", i, rvir[i]);
//		for(j=0; j<r_num; j++)
//		{
//			printf("exc[%d][%d]=%lf\n", i, j ,exc[i][j]);
//		}
//
	}
	for(i=0;i<(m_max-1);i++)
	{
		nc = nc + (Dn_divide_dm[i]*N_cen_average[i] + Dn_divide_dm[i+1]*N_cen_average[i+1])*(M[i+1]-M[i]);
		ns = ns + (Dn_divide_dm[i]*N_sat_average[i] + Dn_divide_dm[i+1]*N_sat_average[i+1])*(M[i+1]-M[i]);
    	ng = ng + (Dn_divide_dm[i]*N_average[i] + Dn_divide_dm[i+1]*N_average[i+1])*(M[i+1]-M[i]);
    	*halo_mass = *halo_mass + (Dn_divide_dm[i]*N_average[i]*M[i]    + Dn_divide_dm[i+1]*N_average[i+1]*M[i+1])   *(M[i+1]-M[i]);
    	*beff      = *beff      + (Dn_divide_dm[i]*N_average[i]*Bias[i] + Dn_divide_dm[i+1]*N_average[i+1]*Bias[i+1])*(M[i+1]-M[i]);
	}
	nc = nc/2.0;
	ns = ns/2.0;
	ng = ng/2.0;
	*fsat = ns/ng;
	*halo_mass = *halo_mass/ng/2.0;
	*beff = *beff/ng/2.0;

//	printf("nc = %lf\n",nc);
//	printf("ns = %lf\n",ns);
//	printf("ng = %lf\n",ng);
//	printf("fsat=%lf\n",fsat);
//	printf("m_mean=%le\n", m_mean);
	
	for(j=0; j<K_num; j++)
	{
		for(i =0; i<(m_max-1); i++)
    	{
    		pcs1[j] = pcs1[j] + (N_cen_average[i]*N_sat_average[i]*Dn_divide_dm[i]*u[i*K_num+j]+\
    			N_cen_average[i+1]*N_sat_average[i+1]*Dn_divide_dm[i+1]*u[(i+1)*K_num+j])*(M[i+1]-M[i]);
    		pss1[j] = pss1[j] + (N_sat_average[i]*N_sat_average[i]*Dn_divide_dm[i]*u[i*K_num+j]*u[i*K_num+j]+\
        		N_sat_average[i+1]*N_sat_average[i+1]*Dn_divide_dm[i+1]*u[(i+1)*K_num+j]*u[(i+1)*K_num+j])*(M[i+1]-M[i]); 
		}
		pcs1[j] = pcs1[j]/ng/ng;
		pss1[j] = pss1[j]/ng/ng/2.0;
		p1[j] =  pcs1[j] + pss1[j];
//		printf("pcs1[%d] = %.8lf\tpss1[%d] = %.8lf\tp1[%d] = %.8lf\n",j, pcs1[j], j, pss1[j], j, p1[j]);
	}
	

	gsl_interp_accel *acc1= gsl_interp_accel_alloc ();
	gsl_interp_accel *acc2= gsl_interp_accel_alloc ();

	gsl_spline *spline1= gsl_spline_alloc (gsl_interp_cspline, K_num);
	gsl_spline *spline2= gsl_spline_alloc (gsl_interp_cspline, K_num);

	gsl_spline_init (spline1, K, Pk, K_num);
	gsl_spline_init (spline2, K, p1, K_num);

	for(mum=0; mum<cut_number; mum++)
	{	
		for(j=r_cut[mum]; j<r_cut[mum+1]; j++)
		{
			for(i=0; i<r_cut_value[mum]; i++)
			{	
				for(num=0; num<node; num++)
				{
					xi_res[num] = Xi[num] + 2*pi*i;
					wi_res[num] = Wi[num];
					corr_mm[j] = corr_mm[j] + wi_res[num]*xi_res[num]/R_corr[j]/R_corr[j]/R_corr[j]*gsl_spline_eval(spline1, xi_res[num]/R_corr[j], acc1)*sin(xi_res[num]);
					corr_1h[j] = corr_1h[j] + wi_res[num]*xi_res[num]/R_corr[j]/R_corr[j]/R_corr[j]*gsl_spline_eval(spline2, xi_res[num]/R_corr[j], acc2)*sin(xi_res[num]);
				}
			}
			corr_mm[j] = corr_mm[j]/2.0/pi/pi;			
			corr_1h[j] = corr_1h[j]/2.0/pi/pi;
		}
	}

	gsl_spline_free (spline1);
	gsl_spline_free (spline2);

	gsl_interp_accel_free (acc1);
	gsl_interp_accel_free (acc2);


	for(j=0; j<r_num; j++)
	{
		for(i=0; i<(m_max-1); i++)
		{
			fexc[j] = fexc[j] + (Dn_divide_dm[i]*N_average[i]*Bias[i]*exc[i*r_num+j] + Dn_divide_dm[i+1]*N_average[i+1]*Bias[i+1]*exc[(i+1)*r_num+j])*(M[i+1]-M[i]);
		}
		fexc[j] = fexc[j]/ng/2.0;
		corr_2hnew[j] = corr_2hnew[j] + fexc[j]*fexc[j]*corr_mm[j];
		corrrrr[j] = corr_1h[j] + corr_2hnew[j];
//		printf("fexc[%d] = %.8lf\n", j, fexc[j]);
//		printf("%d\t%.8lf\t%.8lf\t%.8lf\t%.8lf\t%.8lf\n", j, corr_2h[j], corr_2hold1[j], corr_2hold2[j], corr_2hnew[j], corr_1h[j]+corr_2hnew[j]);
//		printf("%d\t%.8lf\t%.8lf\t%.8lf\n", j, corr_1h[j], corr_2hnew[j], corr_1h[j]+corr_2hnew[j]);
//		printf("%.8lf\n", corr_1h[j]+corr_2hnew[j]);		
	}

//	ft = fopen("corrrrrrrr","r");
//	for(i=0; i<r_test; i++)
//	{
//		fscanf(ft, "%lf", &corrrrr[i]);
//		printf("%.11lf\n", corrrrr[i]);
//	}

	table_pho = gsl_integration_glfixed_table_alloc(node_pho);
	for(num=0; num<number; num++)
	{
		low_limit=dpi_max*num;
		high_limit=dpi_max*(num+1);
		for(i=0; i<node_pho; i++)
		{
			gsl_integration_glfixed_point(low_limit,high_limit,i,&phoxi[num*node_pho+i],&phowi[num*node_pho+i],table_pho);
//			printf("%lf\t%lf\n", phoxi[num*node_pho+i], phowi[num*node_pho+i]);
		}
		
		for(pa = 0; pa<node_pho; pa++)
		{	
			for(pe = 0; pe<Perpendicular_min; pe++)
			{	
				for(i =0; i<node_vel; i++)
				{	
					location=approximation(phoxi[num*node_pho+pa],Perpendicular[pe],Velxi[i]);
					xi_s[num*node_pho+pa][pe]=xi_s[num*node_pho+pa][pe]+Velwi[i]*corrrrr[location]*exp( -pow(Velxi[i],2)/double_sigma_square);

//					printf("%d\n",i);
//					printf("%d\t%d\n",location1,location2);
//					printf("corr1=%.16le\n",corrrrr[location1]);
//					printf("corr2=%.16le\n",corrrrr[location2]);
//					printf("v1=%.16Lf\n",Velocity[i]);
//					printf("v2=%.16Lf\n",Velocity[i+1]);
//					printf("%.16le\n",exp( -pow(Velocity[i],2) /double_sigma_square));
//					printf("%.16le\n",exp( -pow(Velocity[i+1],2) /double_sigma_square));
//					printf("%.16Lf\n",(Velocity[i+1]-Velocity[i]));
				}
				xi_s[num*node_pho+pa][pe]=xi_s[num*node_pho+pa][pe]/(sqrt_2pi*error);	
//				printf("xi_s[%d][%d]=%.14le\n", num*node_pho+pa, pe, xi_s[num*node_pho+pa][pe]);
			}
		}

		for(j=0; j<Perpendicular_min; j++)
		{
			for(i=0; i<node_pho*(num+1); i++)
			{
//				printf("xi_s[%d][%d]=%.14le\n", i, j, xi_s[i][j]);
				w[num][j] = w[num][j] + 2*phowi[i]*xi_s[i][j];
			}
//			printf("%.14lf\n", w[num*perpendicular_max+j]);
		}
	}
	for(i=0; i<number*Perpendicular_min; i++)
	{
		w_analysis[i]=w[i/Perpendicular_min][i%Perpendicular_min];
//		printf("hope%lf\n", w_analysis[i]);
	}
	free(p1);
	free(pcs1);
	free(pss1);
	free(p);
	free(u);
	for(i=0 ;i<node_pho*number; i++)
	{
		free(xi_s[i]);
	}
	for(i=0 ;i<number; i++)
	{
		free(w[i]);
	}	
	return ng;
}







void getMCMC(double *init, int *argc_address, char ***argv_address)
{
	FILE *foutput;
	FILE *fmean, *finvconvar;
	int num,i,j,loop;
	double *ng;
	double *halo_mass;
	double *fsat;	
	double *beff;
	double mag[ndim]={0.01,0.01,3};
	double random0,random1[ndim]={0}, random2[ndim]={0}, rand_delta[ndim]={0};
	double *temp;
	double *delta_w;
	double random_xi;
	double *observe,*w_observe;
	double *covariance_inv;
	double judge;
	double *pos;
	double *w_analysis;
	double *ka_square;
	double *result1,*result2,*result3,*result4,*result5,*result6,*result7,*result8;
	double *result_log_M_min;
	double *result_M1_primer;
	double *result_error;
	double *result_ng;
	double *result_halo_mass;
	double *result_fsat;	
	double *result_beff;	
	double *result_ka;
	int root =0;
	int myid;
		
	temp = 				(double *)malloc(number*Perpendicular_min*sizeof(double));
	observe = 			(double *)malloc(number*perpendicular_max*sizeof(double));
	w_observe = 		(double *)malloc(number*Perpendicular_min*sizeof(double));
	w_analysis =		(double *)malloc(number*Perpendicular_min*sizeof(double));
	covariance_inv =	(double *)malloc(number*Perpendicular_min*number*Perpendicular_min*sizeof(double));
	pos = 				(double *)malloc(walks*ndim*sizeof(double));
	halo_mass = 		(double *)malloc(walks*sizeof(double));
	fsat =				(double *)malloc(walks*sizeof(double));
	beff =				(double *)malloc(walks*sizeof(double));	
	ng = 				(double *)malloc(walks*sizeof(double));
	delta_w = 			(double *)malloc(walks*sizeof(double));
	ka_square = 		(double *)malloc(walks*sizeof(double));
	result1 = 			(double *)malloc(walks*sizeof(double));
	result2 = 			(double *)malloc(walks*sizeof(double));
	result3 = 			(double *)malloc(walks*sizeof(double));
	result4 = 			(double *)malloc(walks*sizeof(double));
	result5 = 			(double *)malloc(walks*sizeof(double));
	result6 = 			(double *)malloc(walks*sizeof(double));
	result7 = 			(double *)malloc(walks*sizeof(double));	
	result8 = 			(double *)malloc(walks*sizeof(double));		
	result_log_M_min =	(double *)malloc(walks*chains*sizeof(double));
	result_M1_primer =	(double *)malloc(walks*chains*sizeof(double));
	result_error =		(double *)malloc(walks*chains*sizeof(double));	
	result_ng = 		(double *)malloc(walks*chains*sizeof(double));
	result_halo_mass =	(double *)malloc(walks*chains*sizeof(double));
	result_fsat = 		(double *)malloc(walks*chains*sizeof(double));
	result_beff = 		(double *)malloc(walks*chains*sizeof(double));			
	result_ka = 		(double *)malloc(walks*chains*sizeof(double));

	for(i=0; i<walks; i++)
	{
		ng[i]       =0;
		halo_mass[i]=0;
		fsat[i]     =0;
		beff[i]     =0; 
		delta_w[i]  =0;
		ka_square[i]=0;
		result1[i]  =0;
		result2[i]  =0;
		result3[i]  =0;
		result4[i]  =0;
		result5[i]  =0;
		result6[i]  =0;
		result7[i]  =0;
		result8[i]  =0;		
	}

	for(i=0; i<walks*chains; i++)
	{
		result_log_M_min[i]=0;
		result_M1_primer[i]=0;
		result_error[i]=0;
		result_ng[i]=0;
		result_halo_mass[i] =0;
		result_fsat[i] =0;	
		result_beff[i] =0;		
		result_ka[i]=0;
	}

	for(i=0; i<number*Perpendicular_min; i++)
	{
		w_analysis[i]=0;
		temp[i]=0;
//		printf("temp[%d]=%lf\n",i, temp[i]);
	}


	fmean = fopen(Wp_mean, "r");
	finvconvar = fopen(Inverse_matrix, "r");
	
	for(i=0; i<number*perpendicular_max; i++)
	{	
		fscanf(fmean, "%lf", &observe[i]);
//		printf("observe[%d]=%lf\n", i, observe[i]);
	}

	for(i=0; i<number*Perpendicular_min; i++)
	{	
		w_observe[i] = observe[(int)(i/Perpendicular_min*perpendicular_max+Start_point-1+i%Perpendicular_min)];
//		printf("%d\n", i/perpendicular_min*perpendicular_max+start_point-1+i%perpendicular_min);
//		printf("w_observe[%d] = %lf\n", i, w_observe[i]);
	}



	for(i=0; i<number*Perpendicular_min; i++)
	{
		for(j=0; j<number*Perpendicular_min; j++)
		{
			fscanf(finvconvar, "%lf", &covariance_inv[i*number*Perpendicular_min+j]);
//			printf("%lf\t", covariance_inv[i][j]);
		}
//		printf("\n");
	}
	fclose(fmean);
	fclose(finvconvar);

	MPI_Init(argc_address, argv_address);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	sleep(myid);
	srand((unsigned) time(NULL));

	for(i=0; i<ndim; i++)
	{
		random0 = rand()/(RAND_MAX+1.0);
		pos[0*ndim+i] = init[i];
//		printf("%lf\t", pos[0][i]);
	}

//	pos[0*ndim+0] = pos[0*ndim+0] + random0*1.0; 
//	pos[0*ndim+1] = pos[0*ndim+1] + random0*1.0;
//	pos[0*ndim+2] = pos[0*ndim+2] + random0*200.0;


	for(loop=0; loop<walks; loop++)
	{
		if(loop!=0)
		{
			for(i=0; i<ndim; i++)
			{	
				random1[i] = rand()/(RAND_MAX+1.0);
				random2[i] = rand()/(RAND_MAX+1.0);
				rand_delta[i] = sqrt(-2.0*pi*log(random1[i]))*cos(2.0*pi*random2[i])*mag[i];
				pos[loop*ndim+i] = pos[(loop-1)*ndim+i] + rand_delta[i];
//				printf("early%lf\t%lf\t%lf\n", rand_delta[i], pos[loop-1][i], pos[loop][i]);
			}
		}

		for(i=0; i<number*Perpendicular_min; i++)
		{
			w_analysis[i]=0;
		}
		ng[loop] = getwp(&pos[loop*ndim], w_analysis, &halo_mass[loop], &fsat[loop], &beff[loop]);
	
		for(i=0; i<number*Perpendicular_min; i++)
		{	
			printf("%.10lf\t", w_analysis[i]);
//			printf("\n");
			for(j=0; j<number*Perpendicular_min; j++)
			{	
//				printf("temp[%d]=%lf\t", i, temp[i]);			
				temp[i] = temp[i]+(w_observe[j]-w_analysis[j])*covariance_inv[j*number*Perpendicular_min+i];	
//				printf("covariance_inv[%d][%d]=%lf\n", j, i,covariance_inv[j][i]);
			}
			delta_w[loop] = delta_w[loop]+temp[i]*(w_observe[i]-w_analysis[i]);
//			printf("temp[%d]=%lf\t", i, temp[i]);
//			printf("w_observe[%d]=%lf\tw_analysis[%d]=%lf\n", i, w_observe[i], i, w_analysis[i]);
//			printf("delta_w[%d]=%lf\n", loop, delta_w[loop]);
		}
		ka_square[loop] = delta_w[loop] + (Number_density-ng[loop])*(Number_density-ng[loop])/(Sigma_ng*Sigma_ng);
//		printf("delta_w[%d]=%lf\n", loop, delta_w[loop]);
//		printf("ka_square[%d]=%lf\n", loop, ka_square[loop]);
			
		for(i=0 ;i<number*Perpendicular_min; i++)
		{
			temp[i] = 0;
			w_analysis[i] = 0; 
		}
		for(i=0; i<ndim; i++)
		{
//			printf("prepos%lf\n",pos[loop][i]);
		}
//		printf("preng%lf\n", ng[loop]);
//		printf("preka%lf\n", ka_square[loop]);
		if(loop!=0)
		{
			judge = exp((ka_square[loop]-ka_square[loop-1])*(-0.5));
//			printf("judge = %le\n", judge);
			if((pos[loop*ndim+0]<10)||(pos[loop*ndim+0]>16)||(pos[loop*ndim+1]<10)||(pos[loop*ndim+1]>16)||\
				(pos[loop*ndim+2]<100)||(pos[loop*ndim+2]>1000))
			{
//				printf("fuck\n");
				for(i=0; i<ndim; i++)
				{
					pos[loop*ndim+i] = pos[(loop-1)*ndim+i];
				}
				ka_square[loop] = ka_square[loop-1];
				ng[loop]        = ng[loop-1];
				halo_mass[loop] = halo_mass[loop-1];
				fsat[loop]      = fsat[loop-1];
				beff[loop]      = beff[loop-1];				
			}
			else if(judge<1)
			{
				random_xi = rand()/(RAND_MAX+1.0);
//				printf("random_xi = %lf\n", random_xi);
				if(random_xi>judge)
				{	
//					printf("fuckagain\n");
					for(i=0; i<ndim; i++)
					{
						pos[loop*ndim+i] = pos[(loop-1)*ndim+i];
					}
					ka_square[loop] = ka_square[loop-1];
					ng[loop]        = ng[loop-1];
					halo_mass[loop] = halo_mass[loop-1];
					fsat[loop]      = fsat[loop-1];	
					beff[loop]      = beff[loop-1];				
				}
			}
		}

	}

	for(loop=0; loop<walks; loop++)
	{	
		result1[loop] = pos[loop*ndim+0];
		result2[loop] = pos[loop*ndim+1];
		result3[loop] = pos[loop*ndim+2];		
		result4[loop] = ng[loop];
		result5[loop] = halo_mass[loop];
		result6[loop] = fsat[loop];
		result7[loop] = beff[loop];		
		result8[loop] = ka_square[loop];
	}

	MPI_Gather(result1, walks, MPI_DOUBLE, result_log_M_min,   walks,   MPI_DOUBLE, root, MPI_COMM_WORLD);
	MPI_Gather(result2, walks, MPI_DOUBLE, result_M1_primer,   walks,   MPI_DOUBLE, root, MPI_COMM_WORLD); 
	MPI_Gather(result3, walks, MPI_DOUBLE, result_error,	   walks,   MPI_DOUBLE, root, MPI_COMM_WORLD); 	
	MPI_Gather(result4, walks, MPI_DOUBLE, result_ng,          walks,   MPI_DOUBLE, root, MPI_COMM_WORLD);
	MPI_Gather(result5, walks, MPI_DOUBLE, result_halo_mass,   walks,   MPI_DOUBLE, root, MPI_COMM_WORLD); 
	MPI_Gather(result6, walks, MPI_DOUBLE, result_fsat,        walks,   MPI_DOUBLE, root, MPI_COMM_WORLD); 
	MPI_Gather(result7, walks, MPI_DOUBLE, result_beff,        walks,   MPI_DOUBLE, root, MPI_COMM_WORLD); 
	MPI_Gather(result8, walks, MPI_DOUBLE, result_ka,          walks,   MPI_DOUBLE, root, MPI_COMM_WORLD); 


	MPI_Barrier(MPI_COMM_WORLD);
	if(myid==0)
	{	
		foutput = fopen("foutput", "w");
		for(i=0 ;i<walks*chains; i++)
		{
			fprintf(foutput, "%.10lf\n", result_log_M_min[i]);
		}
		for(i=0 ;i<walks*chains; i++)
		{
			fprintf(foutput, "%.10lf\n", result_M1_primer[i]);
		}
		for(i=0 ;i<walks*chains; i++)
		{
			fprintf(foutput, "%.10lf\n", result_error[i]);
		}		
		for(i=0 ;i<walks*chains; i++)
		{
			fprintf(foutput, "%.10lf\n", result_ng[i]);
		}
		for(i=0 ;i<walks*chains; i++)
		{
			fprintf(foutput, "%.10lf\n", log10(result_halo_mass[i]));
		}
		for(i=0 ;i<walks*chains; i++)
		{
			fprintf(foutput, "%.10lf\n", result_fsat[i]);
		}
		for(i=0 ;i<walks*chains; i++)
		{
			fprintf(foutput, "%.10lf\n", result_beff[i]);
		}		
		for(i=0 ;i<walks*chains; i++)
		{
			fprintf(foutput, "%.10lf\n", result_ka[i]);
		}
		fclose(foutput);
	}
	MPI_Finalize();

	free(temp);
	free(observe);
	free(w_observe);
	free(w_analysis);
	free(covariance_inv);
	free(pos);
	free(ng);
	free(delta_w);
	free(ka_square);
	free(result1);
	free(result2);
	free(result3);
	free(result4);
	free(result5);
	free(result6);
	free(result7);	
	free(result8);		
	free(result_log_M_min);
	free(result_M1_primer);
	free(result_error);
	free(result_ng);
	free(result_halo_mass);
	free(result_fsat);	
	free(result_beff);		
	free(result_ka);
}

int main(int argc, char **argv)
{
	double omegam;
	double rho;
	struct Linear *li;
	struct Nonlinear *nonli;
	double init[ndim]={1.20249336e+01,   1.34601700e+01,   2.19386074e+02};
	int i, j, num, pa, pe, loop;
	double start,end,cost;
	double xi_res[node]={0},wi_res[node]={0};
	gsl_integration_glfixed_table *table_vel;

	FILE *fl,*fnl;
//	FILE *fv;
	FILE *parafile;
//	fl = fopen("spread_l_noinp_matterpower.dat", "r");
//	fnl = fopen("spread_noinp_matterpower.dat", "r");
	parafile = fopen("infile21","r");
	fscanf(parafile, "%s ", Linear_matter_power_spectrum);	
	fscanf(parafile, "%s",  Nonlinear_matter_power_spectrum);	
	fscanf(parafile, "%d",  &Start_point);
	fscanf(parafile, "%d",  &Perpendicular_min);	
	fscanf(parafile, "%s",  Wp_mean);	
	fscanf(parafile, "%s",  Inverse_matrix);	
	fscanf(parafile, "%lf", &Number_density);	
	fscanf(parafile, "%lf", &Sigma_ng);
	fscanf(parafile, "%lf", &Redshift);
	fscanf(parafile, "%d",  &K_num);	
	fscanf(parafile, "%lf", &Omega_m0);	
	fscanf(parafile, "%lf", &Omega_lambda0);	
/*
	printf("linear_matter_power_spectrum=%s\n",  linear_matter_power_spectrum);	
	printf("nonlinear_matter_power_spectrum=%s\n",  nonlinear_matter_power_spectrum);	
	printf("start_point=%d\n",  start_point);
	printf("perpendicular_min=%d\n",  perpendicular_min);	
	printf("wp_mean=%s\n",  wp_mean);	
	printf("inverse_matrix=%s\n",  inverse_matrix);	
	printf("number_density=%lf\n", number_density);	
	printf("sigma_ng=%lf\n", sigma_ng);	
	printf("z=%lf\n", z);
	printf("k_num=%d\n",  k_num);	
	printf("Omega_m0=%lf\n", Omega_m0);	
	printf("Omega_lambda0=%lf\n", Omega_lambda0);	
*/
	K = (double *)malloc(K_num*sizeof(double));
	Pkl = (double *)malloc(K_num*sizeof(double));
	Pk = (double *)malloc(K_num*sizeof(double));
	Perpendicular = (double *)malloc(Perpendicular_min*sizeof(double));
	li = (struct Linear *)malloc(K_num*sizeof(struct Linear));
	nonli = (struct Nonlinear *)malloc(K_num*sizeof(struct Nonlinear));

	fl = fopen(Linear_matter_power_spectrum, "r");
	if(fl ==NULL)printf("cao");
	fnl = fopen(Nonlinear_matter_power_spectrum, "r");
	if(fnl ==NULL)printf("cao");
	for(i=0; i<K_num; i++)
	{
		fscanf(fl, "%lf%lf", &li[i].kl_r, &li[i].pkl_r);
		K[i] = li[i].kl_r;
		Pkl[i] = li[i].pkl_r;
		fscanf(fnl, "%lf%lf", &nonli[i].k_r, &nonli[i].pk_r);
		Pk[i] = nonli[i].pk_r;
//		printf("%14.4le\t%14.4le\t%14.4le\n", k[i], pkl[i], pk[i]);
	}
	fclose(parafile);
	fclose(fl);
	fclose(fnl);

	for(i=0; i<r_num; i++)
	{
		R_corr[i] = pow(10,(0.01*i -1.0));
//		printf("r_corr[%d]=%lf\n", i, r_corr[i] );
	}

	Rh= Rho_m(0)/(1.9891E30/little_h)*pow((3.08568E22/little_h),3)*Omega_m0;       //(M_sun/h)/(Mpc/h)^3

//	printf("Rho_m(z0) = %.10le\n", Rho_m(z0));
//	printf("%.10le\n", (1.9891E30/little_h));
//	printf("%.10le\n", pow((3.08568E22/little_h),3));
//	printf("rh = %.10le\n", rh);


//	fv = fopen("v_less","r");
/*	for(i=0; i<velocity_num; i++)
	{
//		fscanf(fv, "%lf", &velocity[i]);
		if(i<velocity_num/2.0)
		{
			Velocity[i] = - pow(10, -1 - 0.02*(1+i-velocity_num/2.0) );
		}
		else
		{
			Velocity[i] = pow(10, -1 + 0.02*(i-velocity_num/2.0) );
		}
//		printf("%.11lf\n", velocity[i]);
	}*/

	table_vel = gsl_integration_glfixed_table_alloc(node_vel);
	for(num=0; num<node_vel; num++)
	{
		gsl_integration_glfixed_point(-1500.0,1500.0,num,&Velxi[num],&Velwi[num],table_vel);
//		printf("%lf\t%lf\n", Velxi[num], Velwi[num]);
	}

	for(i=0; i<m_max+1; i++)
	{
		M_add[i] = pow(10,(9.0+0.1*i));
		R_add[i] = pow(M_add[i]/Rh/(4.0/3*pi),(1.0/3));
//		printf("m_add[%d]=%le\tr_add[%d]=%lf\n", i, m_add[i], i, r_add[i]);	
	}
	Sigma_f(Sigma_add,R_add,m_max+1);

//	for(i=0; i<m_max+1; i++)
//	{
//		printf("sigma_add[%d]=%lf\n", i, sigma_add[i]);  
//	}
	
	for(i=0; i<m_max; i++)
	{
		M[i] = pow(10,(9.05+0.1*i));
		R[i] = pow((M[i]/Rh/(4.0/3*pi)),(1.0/3));
//		printf("m[%d]=%le\tr[%d]=%lf\n", i, m[i], i, r[i]);	
	}
	Sigma_f(Sigma,R,m_max);

//	for(i=0; i<m_max; i++)
//	{
//		printf("sigma[%d]=%lf\n", i, sigma[i]);  
//	}                          


	Multiplicity(F);
	bias_f(Bias);
	for(i=0; i<m_max; i++)
	{
    	Dlninvsigma_divide_dm[i] =(log(1.0/Sigma_add[i+1])-log(1.0/Sigma_add[i]))/(M_add[i+1]-M_add[i]);
    	Dn_divide_dm[i] = F[i]*Rh/M[i]*Dlninvsigma_divide_dm[i];
//    	printf("f[%d]=%lf\n",i, f[i]);
//		printf("dlninvsigma_divide_dm[%d] = %le\t dn_divide_dm[%d] = %le\n", i, dlninvsigma_divide_dm[i], i, dn_divide_dm[i]);
//		printf("bias[%d] = %lf\n", i, bias[i]);
	}

	for(i=0; i<parallel_max; i++)
	{
		Parallel[i] = 2.0*i - parallel_max +1;
//		printf("%lf\n", parallel[i]);
	}
	for(i=0; i<Perpendicular_min; i++)
	{
		Perpendicular[i] = pow(10.0,0.1*i+0.1*Start_point-1.1);
//		printf("%lf\n", perpendicular[i]);
	}

//---------------------------------------------------start MCMC----------------------------------------------------------

	start =clock();

	getMCMC(init, &argc, &argv);

	end =clock();
	cost = (end -start)/CLOCKS_PER_SEC ;
	printf("time=%lf\n",cost);

	return 0;
}
