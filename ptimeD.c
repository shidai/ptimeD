// calculate the phase shift between template and simulated (or real) data 
// error of phase can be calculated
// initial guess of phase shift is added
// try to do two-dim template matching
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "ptimeD.h"
#include "T2toolkit.h"
#include "tempo2pred.h"

int main (int argc, char *argv[])
{
	int h,i,j,k;

	/*
	if (argc != 8)
	{
		printf ("Usage: ptime_time -f fname -std tname (-pt tname) -o oname -single (-multi)\n"
	            "Derive the TOAs\n"
	            "fname: data file; tname: templates; oname: output .tim; -std: standard template format; -pt: ptime template;\n"
				"-single: do freq-dependent matching and get one TOA; -multi: do freq-dependent matching and get TOAs for each channel.\n");
	    exit (0);
	}
	*/


	//////////////////////////////////////////////////////
	char fname[128];   // name of data file
	int nstokes;

	int index, n;
	for (i=0;i<argc;i++)
	{
		if (strcmp(argv[i],"-f") == 0)
		{
            index = i + 1;
			n = 0;
			while ( (index + n) < argc && strcmp(argv[index+n],"-std") != 0 && strcmp(argv[index+n],"-pt") != 0 && strcmp(argv[index+n],"-o") != 0 && strcmp(argv[index+n],"-single") != 0)
			{
				n++;
		    }
			//strcpy(fname,argv[++i]);
		}
	}

	// name of different extension of data files
	char name_data[50]; 
	char name_predict[50]; 
	char name_psrparam[50]; 

	char data[] = "[SUBINT]";
	char predict[] = "[T2PREDICT]";
	char psrparam[] = "[PSRPARAM]";

	/////////////////////////////////////////////////////////////////////////////////
	// start to deal with different data file
	for (k = index; k < index + n; k++)
	{
		// get the data file name
		strcpy(fname,argv[k]);
		printf ("%s\n", fname);

		// name of different extension
		strcpy(name_data,fname);
		strcpy(name_predict,fname);
		strcpy(name_psrparam,fname);

		strcat(name_data, data);
		strcat(name_predict, predict);
		strcat(name_psrparam, psrparam);

		////////////////////////////////////////////////////
	
		double psrfreq;
		psrfreq = read_psrfreq(name_psrparam);
		//printf ("psrfreq: %.15lf\n", psrfreq);
	
		double freqRef;
		freqRef = read_obsFreq (name_data);
		//freqRef = 1369.0; // MHz

		double dm;
		dm = readDm(name_data);
		printf ("DM0: %.4lf\n", dm);
	
		////////////////////////////////////////////////
		int nphase;
		int nchn;
		int nsub;
		int npol;
	
		nchn = get_nchan(name_data);	
		npol = get_npol(name_data);	
		nsub = get_subint(name_data);	
		nphase = get_nphase(name_data);	

		//printf ("%d\n", nchn);
		////////////////////////////////////////////////

		////////////////////////////////////////////////////////////////////////////////

		double *p_multi, *p_multi_deDM;
		p_multi = (double *)malloc(sizeof(double)*nchn*npol*nphase);
		p_multi_deDM = (double *)malloc(sizeof(double)*nchn*npol*nphase);

		double *p_temp, *p_temp_deDM;
		p_temp = (double *)malloc(sizeof(double)*npol*nphase);
		p_temp_deDM = (double *)malloc(sizeof(double)*npol*nphase);

		double phaseShift;

		double freq[nchn];
		int n;
		// start to derive toa from different subint
		for (h = 1; h <= nsub; h++)
		{
			// read profiles from data file
			read_prof(name_data,h,p_multi,nphase);
			read_freq(name_data,h,freq,nchn);

			//readfile(argv[2],&n,tt,p_multi);

			// start to derive toas for different channels
			for (i = 0; i < nchn; i++)
			{
				n = 0;
				for (nstokes = 0; nstokes < npol; nstokes++)
				{
					for (j = 0; j < nphase; j++)
					{
						p_temp[n] = p_multi[nstokes*nchn*nphase + i*nphase + j];
						//printf ("%d %lf\n", n, p_temp[n]);
						n++;
					}
				}

				// dedisperse
				phaseShift = phaseShiftDM (dm, freq[i], freqRef, psrfreq);
				deDM (nphase, npol, p_temp, phaseShift, p_temp_deDM);

				n = 0;
				for (nstokes = 0; nstokes < npol; nstokes++)
				{
					for (j = 0; j < nphase; j++)
					{
						p_multi_deDM[nstokes*nchn*nphase + i*nphase + j] = p_temp_deDM[j];
						//printf ("%d %lf\n", j, p_temp_deDM[j]);
						n++;
					}
				}
			}
			write_prof (name_data, h, p_multi_deDM, nphase);
			modify_freq (name_data, h, freqRef, nchn);
		}

		free(p_multi);
		free(p_multi_deDM);
		free(p_temp);
		free(p_temp_deDM);
	}

	return 0;
}
