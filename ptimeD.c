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
	int h,i,j,k,m;

	//////////////////////////////////////////////////////
	char inName[128];   // name of input data file
	char fname[128];   // name of data file
	char ext[128];   // extension of new data file
	char ext0[]="D";   // default extension of new data file
	int nstokes;
	int mode = 0;  // default: creat new file ".D"
	int pmode = 0;  // default: use predictor

	char freqSSBName[128];   // name of freqSSB file
	int freqSSBMode = 0;  // default: no freqSSB read
	char tdisName[128];   // name of tdis file
	int tdisMode = 0;  // default: no tdis read
	int index, n;
	int indexFreqSSB, indexTdis;
	for (i=0;i<argc;i++)
	{
		if (strcmp(argv[i],"-f") == 0)
		{
			index = i + 1;
			n = 0;
			while ( (index + n) < argc && strcmp(argv[index+n],"-e") != 0 && strcmp(argv[index+n],"-np") != 0 && strcmp(argv[index+n],"-tdis") != 0 && strcmp(argv[index+n],"-freqSSB") != 0)
			{
				n++;
			}
			//strcpy(fname,argv[++i]);
		}
		else if (strcmp(argv[i],"-e") == 0)
		{
			strcpy(ext,argv[++i]);
			mode = 1;  // creat new file with new extension
		}
		else if (strcmp(argv[i],"-np") == 0)  // not use predictor
		{
			pmode = 1;  
		}
		else if (strcmp(argv[i],"-freqSSB") == 0)  // not use predictor
		{
			freqSSBMode = 1;
			indexFreqSSB = i + 1;
		}
		else if (strcmp(argv[i],"-tdis") == 0)  // not use predictor
		{
			tdisMode = 1;
			indexTdis = i + 1;
		}
	}

	T2Predictor pred;

	/////////////////////////////////////////////////////////////////////////////////
	// start to deal with different data file
	for (k = 0; k < n; k++)
	{
		// get the data file name
		strcpy(inName,argv[k+index]);
		if (mode == 0)
		{
			createNewfile(inName, fname, ext0);
			printf ("%s\n", fname);
		}
		else
		{
			createNewfile(inName, fname, ext);
			printf ("%s\n", fname);
		}

		////////////////////////////////////////////////////
	
		double freqRef;  // observing central freq at SSB
		double cfreq;  // observing central freq, e.g., 1369 MHz
		freqRef = read_obsFreqSSB (fname);
		cfreq = read_obsFreq (fname);
		//freqRef = 1369.0; // MHz

		double dm;
		dm = readDm(fname);
		printf ("DM0: %.4lf\n", dm);
	
		T2Predictor_Init(&pred);
		  
		int ret;
		if (ret=T2Predictor_ReadFits(&pred,fname))
		{
			printf("Error: unable to read predictor\n");
			exit(1);
		}

		long int imjd, smjd;
		double offs, mjd, subint_offs;
		imjd = stt_imjd(fname);
		smjd = stt_smjd(fname);
		offs = stt_offs(fname);

		double psrfreq;

		////////////////////////////////////////////////
		int nphase;
		int nchn;
		int nsub;
		int npol;
	
		nchn = get_nchan(fname);	
		npol = get_npol(fname);	
		nsub = get_subint(fname);	
		nphase = get_nphase(fname);	

		//printf ("%d\n", nchn);
		////////////////////////////////////////////////

		double freqSSB[nchn];
		double freqRefSSB;
		double tdis1[nchn],tdis2[nchn],shapiro[nchn];
		FILE *fp;

		if (freqSSBMode == 1)
		{
			strcpy(freqSSBName,argv[k+indexFreqSSB]);
			printf ("%s\n", freqSSBName);
  
			if ((fp = fopen(freqSSBName, "r")) == NULL)
		  {
				fprintf (stdout, "Can't open file\n");
				exit(1);
			}

			m = 0;
			while (fscanf(fp, "%lf %lf", &freqSSB[m], &freqRefSSB) == 2)
			{
				//printf ("%lf\n", freqSSB[m]);
				m++;
			}
			  
			if (fclose (fp) != 0)
				fprintf (stderr, "Error closing\n");
		}
		else if (tdisMode == 1)
		{
			strcpy(tdisName,argv[k+indexTdis]);
			printf ("%s\n", tdisName);

			if ((fp = fopen(tdisName, "r")) == NULL)
		  {
				fprintf (stdout, "Can't open file\n");
				exit(1);
			}

			m = 0;
			while (fscanf(fp, "%lf %lf %lf %lf %lf", &freqRefSSB, &tdis1[m],&tdis2[m],&shapiro[m], &freqSSB[m]) == 5)
			{
				//printf ("%lf\n", tdis[m]);
				m++;
			}
			  
			if (fclose (fp) != 0)
				fprintf (stderr, "Error closing\n");
		}

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
			subint_offs = read_offs(fname, h);
			mjd = imjd + (smjd + offs + subint_offs)/86400.0L;
			//printf ("mjd: %lf\n", mjd);
			
			psrfreq = T2Predictor_GetFrequency(&pred,mjd,cfreq);
		
			// read profiles from data file
			read_prof(fname,h,p_multi,nphase);
			read_freq(fname,h,freq,nchn);

			//readfile(argv[2],&n,tt,p_multi);

			// start to derive toas for different channels
			for (i = 0; i < nchn; i++)
			{
				//printf ("Chn%d\n", i);
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
				if (freqSSBMode == 1)
				{
					psrfreq = read_psrfreq (fname);
					phaseShift = phaseShiftDMfreqSSB (freqSSB[i], dm, freqRefSSB, psrfreq);
				}
				else if (tdisMode == 1)
				{
					psrfreq = read_psrfreq (fname);
					phaseShift = phaseShiftDMtdis (tdis1[i], tdis2[i], shapiro[i], psrfreq, freqSSB[i], freqRefSSB, dm);
				}
				else
				{
					//printf ("%lf\n", freq[i]);
					phaseShift = phaseShiftDM (dm, freq[i], pred, mjd, freqRef, psrfreq, pmode);
				}

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
			write_prof (fname, h, p_multi_deDM, nphase);
			//modify_freq (fname, h, freqRef, nchn);
		}

		free(p_multi);
		free(p_multi_deDM);
		free(p_temp);
		free(p_temp_deDM);
	}

	return 0;
}
