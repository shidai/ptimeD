// functions used to derive phase shifts, according to Taylor 1992  
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>
#include "ptimeD.h"
#include "T2toolkit.h"
#include "tempo2pred.h"
//#include "nrutil.h"
#define ITMAX 100000  // Maximum allowed number of iterations.
#define EPS 1.0e-16 // Machine double floating-point precision.
//#define EPS 3.0e-8 // Machine floating-point precision.
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

double A7 (double phase, params param)
//double A7 (double phase, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn)
// calculate function A7 in Taylor 92
//double A7 (int n, double *amp_s, double *amp_p, double *phi_s, double *phi_p, double phase)
{
	double A7=0;
	int i,j;

	for (i = 0; i < param.nchn; i++)
	{
		for (j = 0; j < param.num; j++)
	    {
			A7+=(j+1)*param.a_s[i][j]*param.a_p[i][j]*sin(param.p_s[i][j]-param.p_p[i][j]+(j+1)*phase);
		//printf ("%lf %lf\n", a_s[i], p_s[i]);
		}
	}
	
	return A7;
}

double A7_multi (double phase, params param)
//double A7_multi (double phase, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn, double *rms)
// calculate function A7 in Taylor 92, for multi-frequency channel
//double A7 (int n, double *amp_s, double *amp_p, double *phi_s, double *phi_p, double phase)
{
	double A7=0;
	int i,j;
	double c1, c2, s ;

	for (i = 0; i < param.nchn; i++)
	{
		c1 = 0.0;
		c2 = 0.0;
		s = 0.0;
		for (j = 0; j < param.num; j++)
		{
			c1 += (j+1)*param.a_s[i][j]*param.a_p[i][j]*sin(param.p_s[i][j]-param.p_p[i][j]+(j+1)*phase);
			c2 += param.a_s[i][j]*param.a_p[i][j]*cos(param.p_s[i][j]-param.p_p[i][j]+(j+1)*phase);
			s += param.a_s[i][j]*param.a_s[i][j];
		//printf ("%lf %lf\n", a_s[i], p_s[i]);
		}
		A7 += (c1*c2)/(s*param.rms[i]*param.rms[i]);
	}
	
	return A7;
}

double A9 (double phase, params param)
//double A9 (double phase, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn)
// calculate function A9 in Taylor 92
{
	double A9=0.0, sum=0.0;
	int i,j;

	for (i = 0; i < param.nchn; i++)
	{
	    for (j = 0; j < param.num; j++)
	    {
		    A9+=param.a_s[i][j]*param.a_p[i][j]*cos(param.p_s[i][j]-param.p_p[i][j]+(j+1)*phase);
		    sum+=param.a_s[i][j]*param.a_s[i][j];
		    //printf ("%lf %lf\n", a_s[i], p_s[i]);
		}
	}
	
	A9=A9/sum;

	return A9;
}

int A9_multi (double phase, params param, double *b)
//int A9_multi (double phase, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn, double *b)
{
// calculate function A9 in Taylor 92, for multi-frequency channel
	double A9, sum;
	int i,j;

	for (i = 0; i < param.nchn; i++)
	{
		A9 = 0.0;
		sum = 0.0;
	  for (j = 0; j < param.num; j++)
	  {
		  A9+=param.a_s[i][j]*param.a_p[i][j]*cos(param.p_s[i][j]-param.p_p[i][j]+(j+1)*phase);
		  sum+=param.a_s[i][j]*param.a_s[i][j];
		  //printf ("%lf %lf\n", a_s[i], p_s[i]);
		}
		b[i] = A9/sum;
	}
	
	return 0;
}

int dft_profiles (int N, double *in, fftw_complex *out)
// dft of profiles
{
	//  dft of profiles 
	///////////////////////////////////////////////////////////////////////
	
	//printf ("%lf\n", in[0]);
	//double *in;
	//fftw_complex *out;
	int i;
	double inUse[N];
	fftw_plan p;
	
	//in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
	//out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
	p = fftw_plan_dft_r2c_1d(N, inUse, out, FFTW_MEASURE);

	for (i = 0; i < N; i++)
	{
		inUse[i] = in[i];
	}

	fftw_execute(p);

	fftw_destroy_plan(p);
	//fftw_free(in); 
	//fftw_free(out);
  
	return 0;
}

int error (double phase, double b, double *errphase, double *errb, params param)
//int error (double phase, double b, double *errphase, double *errb, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn)
// calculate the errors of phase, a and b according to Talyor 1992  
{
	double rms,gk,s1,s2;
	int i,j,n;

	gk=0.0;
	s1=0.0;
	s2=0.0;
	n=0;

	for (i = 0; i < param.nchn; i++)
	{
	    for (j = 0; j < param.num; j++)
	    {
		    //gk+=a_p[i]*a_p[i]+b*b*a_s[i]*a_s[i]-2.0*b*a_s[i]*a_p[i]*cos(p_p[i]-p_s[i]-(i+1)*phase)+a*a*1024.0*1024.0-2.0*a*1024.0*a_p[i]*cos(p_p[i])+2.0*a*b*1024.0*a_s[i]*cos(p_s[i]+(i+1)*phase);
		    gk+=param.a_p[i][j]*param.a_p[i][j]+b*b*param.a_s[i][j]*param.a_s[i][j]-2.0*b*param.a_s[i][j]*param.a_p[i][j]*cos(param.p_p[i][j]-param.p_s[i][j]-(j+1)*phase);
		//printf ("%lf %lf\n", a_s[i], p_s[i]);
		s1+=(j+1)*(j+1)*param.a_p[i][j]*param.a_s[i][j]*cos(param.p_p[i][j]-param.p_s[i][j]-(j+1)*phase);
		s2+=param.a_s[i][j]*param.a_s[i][j];
		n++;
		}
	}
	
	rms=sqrt(gk/n);

	(*errphase)=rms/sqrt(2.0*fabs(b*s1));
	(*errb)=rms/sqrt(2.0*fabs(s2));

	return 0;
}

/*
int error_multi (double phase, double *errphase, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn, double *rms, double *bx)
// calculate the errors of phase, a and b according to Talyor 1992  
{
	double s1;
	int i,j,n;

	n=0;
	s1=0.0;

	for (i = 0; i < nchn; i++)
	{
	    for (j = 0; j < num; j++)
	    {
		    s1+=(bx[i]*(j+1)*(j+1)*a_p[i][j]*a_s[i][j]*cos(p_p[i][j]-p_s[i][j]-(j+1)*phase))/(rms[i]*rms[i]);
		    //s2+=(a_s[i][j]*a_s[i][j])/(rms[i]*rms[i]);
		n++;
		}
	}
	
	(*errphase)=1.0/sqrt(2.0*fabs(s1));
	//(*errb)=1.0/sqrt(2.0*s2);

	return 0;
}
*/

int error_multi (double phase, double *errphase, params param)
//int error_multi (double phase, double *errphase, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn, double *rms)
// calculate the errors of phase, a and b according to Talyor 1992  
{
	double s1;
	int i,j;

	s1=0.0;

	double c1, c2, c3, s;

	for (i = 0; i < param.nchn; i++)
	{
		c1 = 0.0;
		c2 = 0.0;
		c3 = 0.0;
		s = 0.0;
		for (j = 0; j < param.num; j++)
		{
			c1 += (j+1)*param.a_s[i][j]*param.a_p[i][j]*sin(param.p_s[i][j]-param.p_p[i][j]+(j+1)*phase);
			c2 += param.a_s[i][j]*param.a_p[i][j]*cos(param.p_s[i][j]-param.p_p[i][j]+(j+1)*phase);
			c3 += (j+1)*(j+1)*param.a_s[i][j]*param.a_p[i][j]*cos(param.p_s[i][j]-param.p_p[i][j]+(j+1)*phase);
			s += param.a_s[i][j]*param.a_s[i][j];
			//s2+=(a_s[i][j]*a_s[i][j])/(rms[i]*rms[i]);
		}
		s1 += (c2*c3-c1*c1)/(s*param.rms[i]*param.rms[i]);
	}
	
	(*errphase)=1.0/sqrt(2.0*fabs(s1));
	//(*errb)=1.0/sqrt(2.0*s2);

	return 0;
}

double zbrent(double (*func)(double phase, params param), double x1, double x2, double tol, params param)
//	Using Brent’s method, find the root of a function func known to lie between x1 and x2. The root, returned as zbrent, will be refined until its accuracy is tol.
{
	int iter;
	double a=x1,b=x2,c=x2,d,e,min1,min2;
	double fa=(*func)(a, param),fb=(*func)(b, param),fc,p,q,r,s,tol1,xm;
	//double fa=(*func)(a, param),fb=(*func)(b, a_s, a_p, p_s, p_p, num, nchn),fc,p,q,r,s,tol1,xm;

	if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
		puts ("Root must be bracketed in zbrent\n");

	fc=fb;
	for (iter=1;iter<=ITMAX;iter++) 
	{
		if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) 
		{
			c=a;   // Rename a, b, c and adjust bounding interval d.
			fc=fa;
			e=d=b-a;
		}
		if (fabs(fc) < fabs(fb)) 
		{
			a=b;
			b=c;
			c=a;
			fa=fb;
			fb=fc;
			fc=fa;
		}

		tol1=2.0*EPS*fabs(b)+0.5*tol;   // Convergence check.
		xm=0.5*(c-b);

		if (fabs(xm) <= tol1 || fb == 0.0) return b;
		if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) 
		{
			s=fb/fa;  // Attempt inverse quadratic interpolation.

			if (a == c) 
			{
				p=2.0*xm*s;
				q=1.0-s;
			} 
			else 
			{
				q=fa/fc;
				r=fb/fc;
				p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
				q=(q-1.0)*(r-1.0)*(s-1.0);
			}
			if (p > 0.0) q = -q;  // Check whether in bounds.

			p=fabs(p);
			min1=3.0*xm*q-fabs(tol1*q);
			min2=fabs(e*q);

			if (2.0*p < (min1 < min2 ? min1 : min2)) 
			{
				e=d;  // Accept interpolation.
				d=p/q;
			} 
			else 
			{
				d=xm; // Interpolation failed, use bisection.
				e=d;
			}
		} 
		else  // Bounds decreasing too slowly, use bisection.
		{
			d=xm;
			e=d;
		}
		a=b;  //  Move last best guess to a.
		fa=fb;
		if (fabs(d) > tol1)     //  Evaluate new trial root.
			b += d;
		else
			b += SIGN(tol1,xm);

		fb=(*func)(b, param);
		//fb=(*func)(b, a_s, a_p, p_s, p_p, num, nchn);
	}

	puts ("Maximum number of iterations exceeded in zbrent\n");

	return 0.0;
}

double zbrent_multi(double (*func)(double phase, params param), double x1, double x2, double tol, params param)
//double zbrent_multi(double (*func)(double phase, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn, double *rms), double x1, double x2, double tol, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn, double *rms)
{
	int iter;
	double a=x1,b=x2,c=x2,d,e,min1,min2;
	double fa=(*func)(a, param),fb=(*func)(b, param),fc,p,q,r,s,tol1,xm;
	//double fa=(*func)(a, a_s, a_p, p_s, p_p, num, nchn, rms),fb=(*func)(b, a_s, a_p, p_s, p_p, num, nchn, rms),fc,p,q,r,s,tol1,xm;

	if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
		puts ("Root must be bracketed in zbrent\n");

	fc=fb;
	for (iter=1;iter<=ITMAX;iter++) 
	{
		if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) 
		{
			c=a;   // Rename a, b, c and adjust bounding interval d.
			fc=fa;
			e=d=b-a;
		}
		if (fabs(fc) < fabs(fb)) 
		{
			a=b;
			b=c;
			c=a;
			fa=fb;
			fb=fc;
			fc=fa;
		}

		tol1=2.0*EPS*fabs(b)+0.5*tol;   // Convergence check.
		xm=0.5*(c-b);

		if (fabs(xm) <= tol1 || fb == 0.0) return b;
		if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) 
		{
			s=fb/fa;  // Attempt inverse quadratic interpolation.

			if (a == c) 
			{
				p=2.0*xm*s;
				q=1.0-s;
			} 
			else 
			{
				q=fa/fc;
				r=fb/fc;
				p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
				q=(q-1.0)*(r-1.0)*(s-1.0);
			}
			if (p > 0.0) q = -q;  // Check whether in bounds.

			p=fabs(p);
			min1=3.0*xm*q-fabs(tol1*q);
			min2=fabs(e*q);

			if (2.0*p < (min1 < min2 ? min1 : min2)) 
			{
				e=d;  // Accept interpolation.
				d=p/q;
			} 
			else 
			{
				d=xm; // Interpolation failed, use bisection.
				e=d;
			}
		} 
		else  // Bounds decreasing too slowly, use bisection.
		{
			d=xm;
			e=d;
		}
		a=b;  //  Move last best guess to a.
		fa=fb;
		if (fabs(d) > tol1)     //  Evaluate new trial root.
			b += d;
		else
			b += SIGN(tol1,xm);

		fb=(*func)(b, param);
		//fb=(*func)(b, a_s, a_p, p_s, p_p, num, nchn, rms);
	}

	puts ("Maximum number of iterations exceeded in zbrent\n");

	return 0.0;
}

/*
int find_peak (int n, double *s, int *position)
{
	int i;
	double temp[n];
	double peak;

	for (i=0;i<n;i++)
	{
		temp[i]=s[i];
	}

	double a,b,c;
	for (i=0;i<n-1;i++)
	{
		a=temp[i];
		b=temp[i+1];
		c=(a>=b ? a : b);

		temp[i+1]=c;
	}
	peak=temp[n-1];

	for (i=0;i<n;i++)
	{
		if (fabs(peak-s[i])<1.0e-3)
		{
			(*position)=i;
		}
	}

	return 0;
}

double find_peak_value (int n, double *s)
{
	int i;
	double temp[n];

	for (i=0;i<n;i++)
	{
		temp[i]=s[i];
	}

	double a,b,c;
	for (i=0;i<n-1;i++)
	{
		a=temp[i];
		b=temp[i+1];
		c=(a>=b ? a : b);

		temp[i+1]=c;
	}

	return temp[n-1];
}
*/

int get_toa (double *s, double *p, double *phasex, double *errphasex, double psrfreq, int nphase, double *rms, double *bx)
// calculate the phase shift between template and simulated (or real) data 
// error of phase can be calculated
// initial guess of phase shift is added
// try to do two-dim template matching
{
    //int nphase=1024;
    int nchn=1;

	// read a std
	
	//puts(argv[1]);
	//puts(argv[2]);
	//double t[nphase*nchn],s[nphase*nchn];
	//int n;

	//readfile(name,&n,t,s);
	//printf ("%d\n", n);
	//puts(argv[1]);

	//////////////////////////////////////////////////////////////////////////
	// simulate data

	//double p[nphase*nchn];
	//double SNR=atof(argv[2]);
	//simulate(n,SNR,s,p);//*/
	
	/////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////
	// dft profile and template
	
	//nchn = n/nphase;
	//printf ("%d\n", nchn);
	int k;  // k=nphase/2

	//double amp_s[nchn][nphase/2],amp_p[nchn][nphase/2];  // elements for calculating A7
	//double phi_s[nchn][nphase/2],phi_p[nchn][nphase/2];
	params param;
	allocateMemoryPtime (&param, nchn, nphase);
	//double amp_s[nchn][NP],amp_p[nchn][NP];  // elements for calculating A7
	//double phi_s[nchn][NP],phi_p[nchn][NP];  // the second dim should be NP, which is large enough for different observations

	preA7(s, p, nphase, nchn, &param);
	//preA7(&k, amp_s, amp_p, phi_s, phi_p, s, p, nphase, nchn);
	//printf ("%d\n", nchn);
	
	// initial guess of the phase
	int peak_s, peak_p;	

	find_peak(0, nphase,s,&peak_s);
	find_peak(0, nphase,p,&peak_p);

	int d, chn;
	double step;
	double ini_phase,up_phase,low_phase;

	d = InitialGuess (s, p, nphase, 1, &chn);
	//d=peak_p-peak_s;
	//printf ("Initial guess: %d\n",d);
	step=2.0*M_PI/(10.0*nphase);

	if (d>=nphase/2)
	{
		if (d>0)
		{
			ini_phase=2.0*M_PI*(nphase-1-d)/nphase;
		}
		else
		{
			ini_phase=-2.0*M_PI*(nphase-1+d)/nphase;
		}
		up_phase=ini_phase+step;
		low_phase=ini_phase-step;
		while (A7(up_phase, param)*A7(low_phase, param)>0.0)
		{
		    up_phase+=step;
		    low_phase-=step;
		}
	}
	else
	{
		ini_phase=-2.0*M_PI*d/nphase;
		up_phase=ini_phase+step;
		low_phase=ini_phase-step;
		while (A7(up_phase, param)*A7(low_phase, param)>0.0)
		{
		    up_phase+=step;
		    low_phase-=step;
		}
	}

  // calculate phase shift, a and b
  double phase,b;
  phase=zbrent(A7, low_phase, up_phase, 1.0e-16, param);
  //phase=zbrent(A7, -1.0, 1.0, 1.0e-16);
  //phase=zbrent(A7, -0.005, 0.005, 1.0e-16);
  b=A9(phase, param);
  //a=A4(b);
	(*bx) = b;

		
	//printf ("Phase shift: %.10lf\n", phase/(2.0*M_PI));
	//printf ("%.10lf %.10lf\n", phase, A7(phase));
	//printf ("%.10lf \n", ((phase/3.1415926)*5.75/2.0)*1.0e+3);  // microseconds
	//printf ("%.10lf \n", b);
	//printf ("%.10lf \n", a);
	//printf ("///////////////////////// \n");
		
	
	// calculate the errors of phase and b
    double errphase, errb;	

	error(phase,b,&errphase,&errb, param);
	//printf ("%.10lf %.10lf\n", ((phase/3.1415926)/(psrfreq*2.0))*1.0e+6, ((errphase/3.1415926)/(psrfreq*2.0))*1.0e+6);  // microseconds
	//printf ("%.10lf %.10lf\n", ((phase/3.1415926)*4.569651/2.0)*1.0e+3, ((errphase/3.1415926)*4.569651/2.0)*1.0e+3);  // microseconds
	//printf ("errphase %.10lf \n", ((errphase/3.1415926)*5.75/2.0)*1.0e+6);
	//printf ("errb %.10lf \n", errb);
	
	(*phasex) = phase;
	(*errphasex) = errphase;

	// calculate the rms
	cal_rms(phase,b,rms, param);

	deallocateMemoryPtime (&param, nchn);
	return 0;
}

int get_toa_multi (char *name_data, char *name_predict, int h, double *s, double *p, double *rms, int nchn, double *phasex, double *errphasex, double psrfreq, int nphase, double *freqout)
//int get_toa_multi (double *s, double *p, double *rms, double *bx, int nchn, double *phasex, double *errphasex, double psrfreq, int nphase)
{
	// get the freq of the subint
	double weight, frequency;
	double freq[nchn], wts[nchn];
	read_freq(name_data, h, freq, nchn);
	read_wts(name_data, h, wts, nchn);
	frequency = 0.0;
	weight = 0.0;

	int z;
	for (z = 0; z < nchn; z++)
	{
		frequency += freq[z]*wts[z];
		weight += wts[z];
	}
	frequency = frequency/weight;
	(*freqout) = frequency;
  //printf ("Frequency is %lf\n", frequency);
	
	int k;  // k=nphase/2

	params param;
	allocateMemoryPtime (&param, nchn, nphase);
	param.rms = rms;

	preA7(s, p, nphase, nchn, &param);
	
	int d, chn;
	double step;
	double ini_phase,up_phase,low_phase;

	d = InitialGuess (s, p, nphase, nchn, &chn);
	//d = InitialGuess (st, pt, nphase, nchn);
	//d=peak_p-peak_s;
	//printf ("Initial guess: %d\n",d);
	step=2.0*M_PI/(10.0*nphase);

	if (d>=nphase/2)
	{
		if (d>0)
		{
			ini_phase=2.0*M_PI*(nphase-1-d)/nphase;
		}
		else
		{
			ini_phase=-2.0*M_PI*(nphase-1+d)/nphase;
		}
		up_phase=ini_phase+step;
		low_phase=ini_phase-step;
		while (A7_multi(up_phase, param)*A7_multi(low_phase, param)>0.0)
		//while (A7_multi(up_phase, amp_s, amp_p, phi_s, phi_p, k, nchn, rms, bx)*A7_multi(low_phase, amp_s, amp_p, phi_s, phi_p, k, nchn, rms,bx)>0.0)
		{
		    up_phase+=step;
		    low_phase-=step;
		}
	}
	else
	{
		ini_phase=-2.0*M_PI*d/nphase;
		up_phase=ini_phase+step;
		low_phase=ini_phase-step;
		while (A7_multi(up_phase, param)*A7_multi(low_phase, param)>0.0)
		//while (A7_multi(up_phase, amp_s, amp_p, phi_s, phi_p, k, nchn, rms, bx)*A7_multi(low_phase, amp_s, amp_p, phi_s, phi_p, k, nchn, rms, bx)>0.0)
		{
		    up_phase+=step;
		    low_phase-=step;
		}
	}
	printf ("initial guess: %lf\n", ini_phase);

  // calculate phase shift, a and b
  double phase;
  phase=zbrent_multi(A7_multi, low_phase, up_phase, 1.0e-16, param);

	double b[nchn];
	A9_multi (phase, param, b);

		
	//printf ("Phase shift: %.10lf\n", phase);
	//printf ("%.10lf \n", ((phase/3.1415926)*5.75/2.0)*1.0e+3);  // microseconds
	//printf ("%.10lf \n", b);
	//printf ("%.10lf \n", a);
	//printf ("///////////////////////// \n");
		
	
	// calculate the errors of phase and b
  double errphase;	

	error_multi(phase, &errphase, param);
	//error_multi(phase, &errphase, amp_s, amp_p, phi_s, phi_p, k, nchn, rms, bx);
	printf ("multi-template\n");
	printf ("Phase shift: %.10lf+-%.10lf\n", ((phase/M_PI)/(psrfreq*2.0))*1.0e+6, ((errphase/M_PI)/(psrfreq*2.0))*1.0e+6);  // microseconds
	//printf ("%.10lf %.10lf\n", ((phase/3.1415926)*4.569651/2.0)*1.0e+3, ((errphase/3.1415926)*4.569651/2.0)*1.0e+3);  // microseconds
	//printf ("errphase %.10lf \n", ((errphase/3.1415926)*5.75/2.0)*1.0e+6);
	//printf ("errb %.10lf \n", errb);
	
	(*phasex) = phase;
	(*errphasex) = errphase;

	deallocateMemoryPtime (&param, nchn);

	return 0;
}

int getToaMultiDM (char *name_data, char *name_predict, int h, double *s, double *p, double *rms, int nchn, double *phasex, double *errphasex, double psrfreq, int nphase, double dm, double *freqout)
//int get_toa_multi (double *s, double *p, double *rms, double *bx, int nchn, double *phasex, double *errphasex, double psrfreq, int nphase)
{
	// get the freq of the subint
	double weight, frequency, obsFreq;
	double freq[nchn], wts[nchn];
	read_freq(name_data, h, freq, nchn);
	read_wts(name_data, h, wts, nchn);
	obsFreq = read_obsFreq (name_data);
	frequency = 0.0;
	weight = 0.0;

	int z;
	for (z = 0; z < nchn; z++)
	{
		frequency += freq[z]*wts[z];
		weight += wts[z];
	}
	frequency = frequency/weight;
	(*freqout) = frequency;
  //printf ("Frequency is %lf\n", frequency);
	
	int k;  // k=nphase/2

	params param;
	allocateMemoryPtime (&param, nchn, nphase);
	param.nfreq = freq;
	param.rms = rms;
	param.psrFreq = psrfreq;
	param.dm = dm;
	param.freqRef = obsFreq;
	//param.freqRef = 1369.0;

	preA7(s, p, nphase, nchn, &param);
	
	int d, chn;
	double step;
	double ini_phase,up_phase,low_phase;

	d = InitialGuess (s, p, nphase, nchn, &chn);
	//printf ("chn: %d\n", chn);

	//int peak_s, peak_p;	
	//find_peak(chn, nphase,s,&peak_s);
	//find_peak(chn, nphase,p,&peak_p);
	//printf ("peak_s: %d\n", peak_s);
	//printf ("peak_p: %d\n", peak_p);
	//d=peak_p-peak_s;

	//printf ("Initial guess: %d\n",d);
	step=2.0*M_PI/(10.0*nphase);

	//printf ("initial guess: %lf\n", (double)(d)/nphase);
	//printf ("delay: %lf\n", nphase*(K*param.dm*param.psrFreq)*(1.0/(param.nfreq[chn]*param.nfreq[chn])-1.0/(param.freqRef*param.freqRef))/(2.0*3.1415926));
	//d = d - (int)(nphase*(K*param.dm*param.psrFreq)*(1.0/(param.nfreq[chn]*param.nfreq[chn])-1.0/(param.freqRef*param.freqRef)));
	//d = d - (int)(d/nphase)*nphase;
	//printf ("Initial guess: %d\n",d);
	if (fabs(d)>=nphase/2)
	{
		if (d>0)
		{
			ini_phase=2.0*M_PI*(nphase-1-d)/nphase;
		}
		else
		{
			ini_phase=-2.0*M_PI*(nphase-1+d)/nphase;
		}
	}
	else
	{
		ini_phase=-2.0*M_PI*d/nphase;
	}

	//printf ("initial guess: %lf\n", ini_phase/(2*M_PI));
  // calculate phase shift, DM, a and b
  double phase, dmFit;
	//printf ("fitDM: Initial guess %f\n",ini_phase);
	miniseNelderMead (&param, ini_phase, &phase, &dmFit);
	//miniseD (&param, ini_phase, &phase, &dmFit);

	// calculate the errors of phase and DM
  double errphase, errDm;	
	covariance (&param, phase, dmFit, &errphase, &errDm);

	// test
	//miniseNelderMeadTest (&param, ini_phase, &phase, &dmFit);
	//errphase = 0.01;
	//errDm = 0.01;

	printf ("multi-template\n");
	//printf ("Phase shift: %.10lf+-%.10lf\n", phase, errphase);  // microseconds
	//printf ("Phase shift: %.10lf+-%.10lf\n", ((phase/3.1415926)/(psrfreq*2.0))*1.0e+6, ((errphase/3.1415926)/(psrfreq*2.0))*1.0e+6);  // microseconds
	printf ("DM: %.10lf   %.10lf\n", dmFit, errDm);
	//printf ("%.10lf %.10lf\n", ((phase/3.1415926)*4.569651/2.0)*1.0e+3, ((errphase/3.1415926)*4.569651/2.0)*1.0e+3);  // microseconds
	//printf ("errphase %.10lf \n", ((errphase/3.1415926)*5.75/2.0)*1.0e+6);
	//printf ("errb %.10lf \n", errb);
	
	(*phasex) = phase;
	(*errphasex) = errphase;

	deallocateMemoryPtime (&param, nchn);

	return 0;
}

int preA7 (double *s, double *p, int nphase, int nchn, params *param)
//int preA7 (int *k, double amp_s[][NP], double amp_p[][NP], double phi_s[][NP], double phi_p[][NP], double *s, double *p, int nphase, int nchn)
// preparation for calculating A7 of Talyor 1992  
{
	// nphase is the dimention of one profile, nchn is number of profiles
	// k is the dimention of amp of one profile 
	param->nchn = nchn;
	int i,j;
	
	/////////////////////////////////////////////////////////////////////////////////

	fftw_complex *out_s;
	fftw_complex *out_p;
	
	out_s = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nphase);
	out_p = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nphase);
	
	double s_temp[nphase];  // store one template and profile
	double p_temp[nphase];  

	int n;
	double r_s[nphase/2],im_s[nphase/2];
	double r_p[nphase/2],im_p[nphase/2];
	for (i = 0; i < nchn; i++)
	{
	    for (j=0;j<nphase;j++)
	    {
		    s_temp[j]=s[i*nphase + j];
		    p_temp[j]=p[i*nphase + j];
	    }

	    dft_profiles(nphase,s_temp,out_s);
	    //printf ("%lf %lf\n", out_s[1][0], out_s[1][1]);

	    dft_profiles(nphase,p_temp,out_p);

	    //double amp_s[N/2],phi_s[N/2];
	    //double amp_p[N/2],phi_p[N/2];

		n = 0;
	    for (j = 0; j <= nphase/2-1; j++)
	    {
		    r_s[j]=out_s[j+1][0];
		    im_s[j]=out_s[j+1][1];
		    r_p[j]=out_p[j+1][0];
		    im_p[j]=out_p[j+1][1];
		    //printf ("%lf %lf\n", r_p[i], im_p[i]);
		    //printf ("%lf %lf\n", out_s[i][0], out_s[i][1]);
		    n++;
	    }
	    //printf ("%d\n", n);
	    //printf ("%d %d\n", nphase, nchn);

	    for (j = 0; j < n; j++)
	    {
		    param->a_s[i][j]=sqrt(r_s[j]*r_s[j]+im_s[j]*im_s[j]);
		    param->a_p[i][j]=sqrt(r_p[j]*r_p[j]+im_p[j]*im_p[j]);
		    param->p_s[i][j]=atan2(im_s[j],r_s[j]);
		    param->p_p[i][j]=atan2(im_p[j],r_p[j]);
		    //printf ("%lf %lf %lf\n", r_s[i], im_s[i], amp_s[i]);
		    //printf ("%lf %lf %lf\n", r_p[i], im_p[i], amp_p[i]);
		    //printf ("%lf\n", amp_s[i]);
		    //printf ("%lf\n", amp_p[i]);
	    }
	}
	//(*k)=n;
	param->num = n;

	fftw_free(out_s); 
	fftw_free(out_p); 

	return 0;
}

int cal_rms (double phase, double b, double *rms, params param)
// calculate the rms of each subchannel  
{
	double gk;
	int i,j,n;

	gk=0.0;
	n=0;

	for (i = 0; i < param.nchn; i++)
	{
	    for (j = 0; j < param.num; j++)
	    {
		    //gk+=a_p[i]*a_p[i]+b*b*a_s[i]*a_s[i]-2.0*b*a_s[i]*a_p[i]*cos(p_p[i]-p_s[i]-(i+1)*phase)+a*a*1024.0*1024.0-2.0*a*1024.0*a_p[i]*cos(p_p[i])+2.0*a*b*1024.0*a_s[i]*cos(p_s[i]+(i+1)*phase);
		    gk+=param.a_p[i][j]*param.a_p[i][j]+b*b*param.a_s[i][j]*param.a_s[i][j]-2.0*b*param.a_s[i][j]*param.a_p[i][j]*cos(param.p_p[i][j]-param.p_s[i][j]-(j+1)*phase);
		//printf ("%lf %lf\n", a_s[i], p_s[i]);
		n++;
		}
	}
	
	(*rms)=sqrt(gk/n);

	return 0;
}

/*
int simulate (int n, double SNR, double *s, double *p)
// simulate pulse profiles, adding white noise; return simulated profiles
{
	// simulate a profile with white noise
	///////////////////////////////////////////////////////////////////////
	// initialize gsl 
	
	int i;
	const gsl_rng_type * T;
	gsl_rng * r;

	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc (T);
  
	////////////////////////////////////////////////////////////////////////
	//  determine the amplitude of white noise according to SNR
	
	double scale;   // the scale multiply to white noise to get certain SNR
	double amp_noise, noise[n];
	double peak_s;

	for (i=0;i<n;i++)
	{
		noise[i]=gsl_ran_gaussian(r,1.0);
		amp_noise+=noise[i]*noise[i];
	}
	
	amp_noise=sqrt(amp_noise/n);

	peak_s=find_peak_value(n,s);   // find the peak flux of the std
    //printf ("peak of std: %g\n", peak_s);

	scale=peak_s/(SNR*amp_noise);
    //printf ("%g\n", scale);
	
	//////////////////////////////////////////////////////////////////////////
	//  add noise to std ==> p

	for (i=0;i<n;i++)
	{
		p[i]=(s[i]+scale*noise[i]);
	}

	
	double peak_p;
	peak_p=find_peak(n,p);  // find the peak flux of the profile
    //printf ("peak of profile: %g\n", peak_p);

	//  normalize the std and profile
	
	for (i=0;i<n;i++)
	{
		p[i]=p[i]/peak_p;
		s[i]=s[i]/peak_s;
		//printf ("%g %g\n", s[i], p[i]);
	}

	//double rms=0.0;
	//int m=0;

	//for (i=300;i<700;i++)
	//{
    //	rms+=p[i]*p[i];
    //	m++;
	//}

	//rms=sqrt(rms/m);
	//printf("rms is: %f\n", rms);

	gsl_rng_free (r);
  
	return 0;
}
*/

int form_toa_multi (char *name_data, char *name_predict, int subint, int nchn, long int imjd, long int smjd, double offs, double phase, double e_phase, long double *t, long double *e_dt, double frequency)
{
	int h;
	h = subint;

	long double dt;  
	double offset;   // offset of each subint
	long double mjd0;  // the mjd of each subint
	T2Predictor pred;
	int ret;
	double period; 

	// transform phase shift to TOAs

	// get the period
    print_t2pred(name_predict);   // output t2pred.dat

	T2Predictor_Init(&pred);  // prepare the predictor
	if (ret=T2Predictor_Read(&pred,(char *)"t2pred.dat"))
    {
		printf("Error: unable to read predictor\n");
		exit(1);
	}

	// get the offset of each subint
	offset = read_offs(name_data, h);

	// get the period at mjd0
	mjd0 = (long double)(imjd) + ((long double)(smjd) + (long double)(offs) + (long double)(offset))/86400.0L;
	//printf ("imjd is %ld \n", imjd);
	//printf ("mjd0 is %.15Lf \n", mjd0);

	period = 1.0/T2Predictor_GetFrequency(&pred,mjd0,frequency);
    //printf ("Period is %.15lf\n", period);
	
	// transform phase shift to time shift
    //dt = (phase/M_PI)*period/2.0;
    //e_dt = (e_phase/M_PI)*period/2.0;
    dt = ((long double)(phase)/M_PI)*((long double)(period))/2.0L;
    (*e_dt) = ((long double)(e_phase)/M_PI)*((long double)(period))/2.0L;
    //printf ("dt is %.10Lf +/- %.10Lf\n", dt, e_dt);

	// calculate the TOA
    (*t) = (long double)(imjd) + ((long double)(smjd) + (long double)(offs) - (long double)(dt) + (long double)(offset))/86400.0L;
    //t = imjd;
		
    //printf ("offset is %lf\n", offset);
	//fprintf (fp, "%s  %lf  %.15Lf  %Lf  7\n", fname, frequency, t, e_dt*1e+6);

	return 0;
}

int form_toa (char *name_data, char *name_predict, int subint, int chn, int nchn, long int imjd, long int smjd, double offs, double phase, double e_phase, long double *t, long double *e_dt, double *freqout)
// chn is the channel to form toa
// nchn is the total number of subchn
{
	int h;
	h = subint;

	long double dt;  
	double offset;   // offset of each subint
	long double mjd0;  // the mjd of each subint
	T2Predictor pred;
	int ret;
	double period, frequency;
	double freq[nchn];

	// transform phase shift to TOAs
	// get the freq of the subint
	read_freq(name_data, h, freq, nchn);

	frequency = freq[chn];
	(*freqout) = frequency;
    printf ("Frequency is %lf\n", frequency);

	// get the period
    print_t2pred(name_predict);   // output t2pred.dat

	T2Predictor_Init(&pred);  // prepare the predictor
	if (ret=T2Predictor_Read(&pred,(char *)"t2pred.dat"))
    {
		printf("Error: unable to read predictor\n");
		exit(1);
	}

	// get the offset of each subint
	offset = read_offs(name_data, h);

	// get the period at mjd0
	mjd0 = (long double)(imjd) + ((long double)(smjd) + (long double)(offs) + (long double)(offset))/86400.0L;
	printf ("imjd is %ld \n", imjd);
	printf ("mjd0 is %.15Lf \n", mjd0);

	period = 1.0/T2Predictor_GetFrequency(&pred,mjd0,frequency);
    printf ("Period is %.15lf\n", period);
	
	// transform phase shift to time shift
    //dt = (phase/M_PI)*period/2.0;
    //e_dt = (e_phase/M_PI)*period/2.0;
    dt = ((long double)(phase)/M_PI)*((long double)(period))/2.0L;
    (*e_dt) = ((long double)(e_phase)/M_PI)*((long double)(period))/2.0L;
    printf ("dt is %.10Lf +/- %.10Lf\n", dt, e_dt);

	// calculate the TOA
    (*t) = (long double)(imjd) + ((long double)(smjd) + (long double)(offs) - (long double)(dt) + (long double)(offset))/86400.0L;
    //t = imjd;
		
    printf ("offset is %lf\n", offset);
	//fprintf (fp, "%s  %lf  %.15Lf  %Lf  7\n", fname, frequency, t, e_dt*1e+6);

	return 0;
}

int find_peak (int n0, int n, double *s, int *position)
{
	int i;
	double temp[n];
	double peak;

	for (i = n*n0; i < n*n0+n; i++)
	{
		temp[i-n*n0] = s[i];
	}

	double a, b, c;
	for (i = 0; i < n-1; i++)
	{
		//a = fabs(temp[i]);
		//b = fabs(temp[i+1]);
		a = temp[i];
		b = temp[i+1];
		c = (a >= b ? a : b);

		temp[i+1] = c;
	}
	peak = temp[n-1];

	for (i = n*n0; i < n*n0+n; i++)
	{
		if (fabs(peak-s[i]) < 1.0e-3)
		{
			(*position) = i;
		}
	}

	return 0;
}

double find_peak_value (int n, double *s)
{
	int i;
	double temp[n];

	for (i = 0; i < n; i++)
	{
		temp[i] = s[i];
	}

	double a, b, c;
	for (i = 0; i < n-1; i++)
	{
		a = temp[i];
		b = temp[i+1];
		//a = fabs(temp[i]);
		//b = fabs(temp[i+1]);
		//c = (fabs(a) >= fabs(b) ? a : b);
		c = (a >= b ? a : b);

		temp[i+1] = c;
	}

	return temp[n-1];
}

int corr (double *s, double *p, int nphase)
{
	/*
	FILE *fp1, *fp2;

	if ((fp1 = fopen(argv[1], "r")) == NULL)
	{
		fprintf (stdout, "Can't open file\n");
		exit(1);
	}

	if ((fp2 = fopen(argv[2], "r")) == NULL)
	{
		fprintf (stdout, "Can't open file\n");
		exit(1);
	}

	float x1[1024], x2[1024];
	int i = 0;
	while (fscanf (fp1, "%f", &x1[i]) == 1)
	{
		i++;
	}

	i = 0;
	while (fscanf (fp2, "%f", &x2[i]) == 1)
	{
		i++;
	}
	*/

	double r[nphase];
	int i, j;
	for (j = 0; j < nphase; j++)
	{
		r[j] = 0.0;
		for (i = 0; i < nphase; i++)
		{
			if ((i+j) > (nphase-1))
			{
				r[j] += p[i]*s[i+j-(nphase-1)];
			}
			else
			{
				r[j] += p[i]*s[i+j];
			}
			//printf ("%f %f\n", x1[i], x2[i]);
		}
	}

	int shift;
	find_peak (0, nphase, r,  &shift);
	/*
	for (j = 0; j < 1024; j++)
	{
		printf ("%f\n", r[j]);
	}
	*/

	return -shift;
}

int def_off_pulse (int nphase, double *in, double frac_off)
// define the off pulse region based on I, return the starting index of off pulse region
// using frac_off to calculate the off pulse region
{
	int n = nphase;
	int num_off = (int)(n*frac_off);
	int i,j;
	double small;
	double temp;
	int index = 0;

	for (i = 0; i < n; i++)
	{
		if (i == 0)
		{
			small = 0.0;
			for(j = 0; j < num_off; j++)
			{
				small += (in[j]+30000.0)*(in[j]+30000.0);  // make all numbers positive
			}
			small = sqrt(small/num_off);
		}
			
		temp = 0.0;
		for(j = 0; j < num_off; j++)
		{
			if ((i+j) > n-1)
			{
				temp += (in[(i+j)-(n-1)]+30000.0)*(in[(i+j)-(n-1)]+30000.0);
			}
			else 
			{
				temp += (in[i+j]+30000.0)*(in[i+j]+30000.0);
			}
		}
		temp = sqrt(temp/num_off);

		small = (temp <= small ? temp : small);
		index = (temp <= small ? i : index);
		//printf ("%d %lf %lf\n", index, small, ave);
	}

	return index;
}

int off_pulse (int nphase, int index, double *in, double *out, double frac_off)
// get the off_pulse region
{
	int n = nphase;
	int num_off = (int)(n*frac_off);
	int i;

	for (i = 0; i < num_off; i++)
	{
		if ((index+i) > n-1)
		{
			out[i] = in[(index+i)-(n-1)];
		}
		else 
		{
			out[i] = in[index+i];
		}
	}

	return 0;
}

int remove_baseline (double *in, int index, double frac_off, int n, double *out)
{
	// define the off_pulse range, frac_off is the fraction of the phase
	// index is the starting point of the off_pulse range
	int num_off = (int)(n*frac_off);
	double off_0[num_off];

	off_pulse (n, index, in, off_0, frac_off);

	int i;
	double baseline = 0.0;
    for (i = 0; i < num_off; i++)
    {
        baseline += off_0[i];
        //average_s += s_off[i];
    }
	baseline = baseline/num_off;

    //printf ("the baseline of std is: %lf \n", baseline);
    //printf ("average is: %lf %lf\n", average, average_s);

	for (i = 0; i < n; i++)
	{
		out[i] = (in[i]-baseline);
		//s_norm[i] = (s[i]-baseline)/(s_peak-baseline);
	}
	
	return 0;
}

int pre_diff (double *s, int nphase, int index, double frac_off, double *s_out)
{
	int n = nphase;
	
	// remove the baseline
	remove_baseline (s, index, frac_off, n, s_out);

    return 0;
}


int InitialGuess (double *s, double *p, int nphase, int nchn, int *chn) 
{
	int index;
	double frac_off = 0.05;  // set to be 0.05

	double ptemp[nphase];
	double temp[nchn];
	int i, h;
	for (i = 0; i < nchn; i++)
	{
		for (h = 0; h < nphase; h++)
		{
			ptemp[h] = p[i*nphase+h];
		}
		int x;
		x = def_off_pulse (nphase, ptemp, frac_off);

		double ptemp_out[nphase];
		pre_diff (ptemp, nphase, x, frac_off, ptemp_out);

		temp[i] = find_peak_value(nphase,ptemp_out);
	}

	int peak;
	find_peak (0, nchn, temp, &peak);
	//printf ("%d\n",peak);
	(*chn) = peak;

	double p_use[nphase];
	double s_use[nphase];
	for (h = 0; h < nphase; h++)
	{
		p_use[h] = p[peak*nphase+h];
		s_use[h] = s[peak*nphase+h];
	}

	// remove the baseline of template
	index = def_off_pulse (nphase, s_use, frac_off);

	double s_out[nphase];
	pre_diff (s_use, nphase, index, frac_off, s_out);

	// remove the baseline of profile
	index = def_off_pulse (nphase, p_use, frac_off);

	double p_out[nphase];
	pre_diff (p_use, nphase, index, frac_off, p_out);

	// Guess the phase shift
	int d;
	d = corr (s_out, p_out, nphase);

	return d;
}

int preA7_QUV (double *p, int nphase, double *real_p, double *ima_p)
// preparation for calculating A7 of Talyor 1992  
{
	// nphase is the dimention of one profile, nchn is number of profiles
	// k is the dimention of amp of one profile 
	int i,j;
	
	/////////////////////////////////////////////////////////////////////////////////

	fftw_complex *out_p;
	
	out_p = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nphase);
	
	double p_temp[nphase];  // store one template and profile

	for (j=0;j<nphase;j++)
	{
	  p_temp[j]=p[j];
	}

	dft_profiles(nphase,p_temp,out_p);

	//double amp_s[N/2],phi_s[N/2];
	//double amp_p[N/2],phi_p[N/2];

	for (j = 0; j < nphase/2+1; j++)                                                  
	{                                                                      
		real_p[j]=out_p[j][0];                                             
		ima_p[j]=out_p[j][1];                                              
	}

	fftw_free(out_p); 

	return 0;
}

int rotate (int N, double *real_p, double *real_p_rotate, double *ima_p, double *ima_p_rotate, double rot)
{
	// k is the dimention of amp, N is the dimention of s
	int i;

	// for substraction 
	double amp,cosina,sina;
	for (i=0;i<N/2+1;i++)
	{
		// calculate the sin(phi) and cos(phi) of the profile
		amp=sqrt(real_p[i]*real_p[i]+ima_p[i]*ima_p[i]);
		cosina=real_p[i]/amp;
		sina=ima_p[i]/amp;

		// rotate profile
		real_p_rotate[i]=amp*(cosina*cos(-i*rot)-sina*sin(-i*rot));
		ima_p_rotate[i]=amp*(sina*cos(-i*rot)+cosina*sin(-i*rot));
		//real_p_rotate[i]=amp*(cosina*cos(-i*rot*M_PI)-sina*sin(-i*rot*M_PI));
		//ima_p_rotate[i]=amp*(sina*cos(-i*rot*M_PI)+cosina*sin(-i*rot*M_PI));
		//real_p_rotate[i]=amp*(cosina*cos(-i*M_PI)-sina*sin(-i*M_PI));
		//ima_p_rotate[i]=amp*(sina*cos(-i*M_PI)+cosina*sin(-i*M_PI));
		
	}

	return 0;
}

int align (int N, double phase, double b, double a, double *real_p, double *real_p_align, double *ima_p, double *ima_p_align, double rotate)
{
	// k is the dimention of amp, N is the dimention of s
	int i;

	// for substraction 
	double amp,cosina,sina;
	for (i=0;i<N/2+1;i++)
	{
		// calculate the sin(phi) and cos(phi) of the profile
		amp=sqrt(real_p[i]*real_p[i]+ima_p[i]*ima_p[i]);
		cosina=real_p[i]/amp;
		sina=ima_p[i]/amp;

		// add phase shift to the profile, phase
		//real_p_align[i]=amp*(cosina)/b;
		//ima_p_align[i]=amp*(sina)/b;
		//real_p_align[i]=amp*(cosina*cos(-i*phase)-sina*sin(-i*phase));
		//ima_p_align[i]=amp*(sina*cos(-i*phase)+cosina*sin(-i*phase));
		//real_p_align[i]=amp*(cosina*cos(-i*phase)-sina*sin(-i*phase))/b;
		//ima_p_align[i]=amp*(sina*cos(-i*phase)+cosina*sin(-i*phase))/b;
		real_p_align[i]=(amp*(cosina*cos(-i*(phase+rotate*M_PI))-sina*sin(-i*(phase+rotate*M_PI))))/b;
		ima_p_align[i]=(amp*(sina*cos(-i*(phase+rotate*M_PI))+cosina*sin(-i*(phase+rotate*M_PI))))/b;
		
	}

	return 0;
}

int inverse_dft (double *real_p, double *ima_p, int ncount, double *p_new)
{
	double *dp;
    fftw_plan plan;
	fftw_complex *cp;

    dp = (double *)malloc(sizeof (double) * ncount);
	cp = (fftw_complex *)fftw_malloc(sizeof (fftw_complex) * ncount);
	memset(dp, 0, sizeof (double) * ncount);
	memset(cp, 0, sizeof (fftw_complex) * ncount);

	// initialize the dft...
	double *dp_t;
    fftw_plan plan_t;
	fftw_complex *cp_t;

    dp_t = (double *)malloc(sizeof (double) * ncount);
	cp_t = (fftw_complex *)fftw_malloc(sizeof (fftw_complex) * ncount);
	memset(dp_t, 0, sizeof (double) * ncount);
	memset(cp_t, 0, sizeof (fftw_complex) * ncount);

	int i;
    double real,ima,amp,cosina,sina;

	for (i = 0; i < ncount; i++)
	{
		if (i < ncount/2+1)
		{
            real = real_p[i];
            ima = ima_p[i];
			amp = sqrt(real*real+ima*ima);
			cosina = real/amp;
			sina = ima/amp;

			cp[i][0] = amp*(cosina);
			cp[i][1] = amp*(sina);
			//cp[i][0] = amp*(cosina*cos(-i*3.1415926)-sina*sin(-i*3.1415926));
			//cp[i][1] = amp*(sina*cos(-i*3.1415926)+cosina*sin(-i*3.1415926));
			//cp[i][0]=real_s[i]-real_p[i];
			//cp[i][1]=ima_s[i]-ima_p[i];
			//cp[i][0]=-real_s[i]+real_p[i];
			//cp[i][1]=-ima_s[i]+ima_p[i];
			cp_t[i][0] = real_p[i];
			cp_t[i][1] = ima_p[i];
			//cp[i][0]=real_p[i];
			//cp[i][1]=ima_p[i];
		}
		else
		{
			cp[i][0]=0.0;
			cp[i][1]=0.0;
			cp_t[i][0]=0.0;
			cp_t[i][1]=0.0;
		}
	}

    plan_t = fftw_plan_dft_c2r_1d(ncount, cp_t, dp_t, FFTW_MEASURE);

    fftw_execute(plan_t);

    fftw_destroy_plan(plan_t);

	/////////////////////////////////////////////////////////////////

    plan = fftw_plan_dft_c2r_1d(ncount, cp, dp, FFTW_MEASURE);

    fftw_execute(plan);

    fftw_destroy_plan(plan);

	for (i = 0; i < ncount; i++)
	{
		p_new[i] = dp[i]/ncount;  // normalized by the ncount
		//printf ("%lf\n", p_new[i]);
	}

	return 0;
}

int allocateMemoryPtime (params *param, int nchn, int nphase)
{
	param->a_s = (double **)malloc(sizeof(double *)*nchn);
	param->a_p = (double **)malloc(sizeof(double *)*nchn);
	param->p_s = (double **)malloc(sizeof(double *)*nchn);
	param->p_p = (double **)malloc(sizeof(double *)*nchn);
	//printf ("test\n");

	int i;
	for (i = 0; i < nchn; i++)
	{
		param->a_s[i] = (double *)malloc(sizeof(double)*nphase);
		param->a_p[i] = (double *)malloc(sizeof(double)*nphase);
		param->p_s[i] = (double *)malloc(sizeof(double)*nphase);
		param->p_p[i] = (double *)malloc(sizeof(double)*nphase);
	}
}

int deallocateMemoryPtime (params *param, int nchn)
{
	int i;
	for (i = 0; i < nchn; i++)
	{
		free(param->a_s[i]);
		free(param->a_p[i]);
		free(param->p_s[i]);
		free(param->p_p[i]);
	}
		
	free(param->a_s);
	free(param->a_p);
	free(param->p_s);
	free(param->p_p);
}

int deDM (int nphase, int npol, double *in, double phaseShift, double *out)
// de-disperse 
{
	//printf ("npol: %d\n", npol);
	int i, j;
	
	double I_in[nphase], Q_in[nphase], U_in[nphase], V_in[nphase];
	double I_out[nphase], Q_out[nphase], U_out[nphase], V_out[nphase];
	
	for (i = 0; i < npol; i++)
	{
		for (j = 0; j < nphase; j++)
		{
			if (i == 0)
			{
				I_in[j] = in[i*nphase+j];
			}
			else if (i == 1)
			{
				Q_in[j] = in[i*nphase+j];
			}
			else if (i == 2)
			{
				U_in[j] = in[i*nphase+j];
			}
			else if (i == 3)
			{
				V_in[j] = in[i*nphase+j];
			}
			else
			{
				printf ("Wrong npol!\n");
			}
		}
	}

	double I_in_real[nphase/2+1], I_in_ima[nphase/2+1];
	double Q_in_real[nphase/2+1], Q_in_ima[nphase/2+1];
	double U_in_real[nphase/2+1], U_in_ima[nphase/2+1];
	double V_in_real[nphase/2+1], V_in_ima[nphase/2+1];

	preA7_QUV (I_in, nphase, I_in_real, I_in_ima);

	if (npol == 4)
	{
		preA7_QUV (Q_in, nphase, Q_in_real, Q_in_ima);
  	preA7_QUV (U_in, nphase, U_in_real, U_in_ima);
  	preA7_QUV (V_in, nphase, V_in_real, V_in_ima);
	}

	double I_out_real[nphase/2+1], I_out_ima[nphase/2+1];
	double Q_out_real[nphase/2+1], Q_out_ima[nphase/2+1];
	double U_out_real[nphase/2+1], U_out_ima[nphase/2+1];
	double V_out_real[nphase/2+1], V_out_ima[nphase/2+1];
	rotate (nphase, I_in_real, I_out_real, I_in_ima, I_out_ima, phaseShift);
	if (npol == 4)
	{
		rotate (nphase, Q_in_real, Q_out_real, Q_in_ima, Q_out_ima, phaseShift);
		rotate (nphase, U_in_real, U_out_real, U_in_ima, U_out_ima, phaseShift);
		rotate (nphase, V_in_real, V_out_real, V_in_ima, V_out_ima, phaseShift);
	}

	inverse_dft (I_out_real, I_out_ima, nphase, I_out);
	if (npol == 4)
	{
		inverse_dft (Q_out_real, Q_out_ima, nphase, Q_out);
		inverse_dft (U_out_real, U_out_ima, nphase, U_out);
		inverse_dft (V_out_real, V_out_ima, nphase, V_out);
	}

	for (i = 0; i < npol; i++)
	{
		for (j = 0; j < nphase; j++)
		{
			if (i == 0)
			{
				out[i*nphase+j] = I_out[j];
				//printf ("%d %lf %lf\n", j, I_in[j], I_out[j]);
			}
			else if (i == 1)
			{
				out[i*nphase+j] = Q_out[j];
			}
			else if (i == 2)
			{
				out[i*nphase+j] = U_out[j];
			}
			else if (i == 3)
			{
				out[i*nphase+j] = V_out[j];
			}
			else
			{
				printf ("Wrong npol!\n");
			}
		}
	}

	return 0;
}

