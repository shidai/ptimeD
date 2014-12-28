// read the predictor from PSRFITs files
#include "tempo2pred.h"
#include "tempo2pred_int.h"
#include "ptimeD.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifdef sun
#include <sunmath.h>
#endif

#ifndef M_PIl
#define M_PIl 3.14159265358979323846264338327950288L
#endif

int T2Predictor_ReadFits(T2Predictor *t2p, char *fname)
{
  fitsfile *f;
	int status; 
	status = 0;

	if ( fits_open_file(&f, fname, READONLY, &status) )          // open the file
	{
		printf( "error while openning file\n" );
	}

	// move to the Predictor
	fits_movnam_hdu(f,BINARY_TBL,(char *)"T2PREDICT",0,&status);

  int ret;
  //printf("step1\n");
  if (!f)
    return -1;
  //printf("Got here\n");
  ret = T2Predictor_FReadFits(t2p, f);
  //printf("and here %d\n",ret);
	
	if ( fits_close_file(f, &status) )
	{
		printf( " error while closing the file\n " );
	}

  return ret;
}
  
int T2Predictor_FReadFits(T2Predictor *t2p, fitsfile *f)
{
  // determine the kind of file we're dealing with by trial and error
  if (ChebyModelSet_ReadFits(&t2p->modelset.cheby, f)==0)
  {
    t2p->kind = Cheby;
    return 0;
  }
 
  //if (T1PolycoSet_Read(&t2p->modelset.t1, f)==0)
  //{
  //  t2p->kind = T1;
  //  return 0;
  //}

	return -1;
}
  
int ChebyModelSet_ReadFits(ChebyModelSet *cms, fitsfile *f)
{
	char **line;
	line = (char **)malloc(sizeof(char *));
	line[0] = (char *)malloc(sizeof(char)*1024);
  char keyword[64];

	int colnum;
	int status = 0;
	char nval[]="NULL";
	int anynull = 0;

	fits_get_colnum(f,CASEINSEN,"PREDICT",&colnum,&status);

	fits_read_col(f, TSTRING, colnum, 1, 1, 1, nval, line, &anynull, &status);           // read the column

	//int nchar = strlen(line[0]);

  int iseg;
  int ret;
  //printf("Got here b\n");
  if (sscanf(line[0], "%s %d", keyword, &cms->nsegments)!=2)
    return -1;
  //printf("Got here c\n");
  if (strcasecmp(keyword, "ChebyModelSet"))
    return -1;
  //printf("Got here d\n");
  if (!(cms->segments = (ChebyModel *)malloc(cms->nsegments*sizeof(ChebyModel))))
    {
      printf("Unable to allocate memory in reading predictor\n");
      exit(1);
    }
  //printf("Got here e %d\n",cms->nsegments);
	int row0 = 2;
  for (iseg=0; iseg < cms->nsegments ; iseg++)
	{
		//printf ("%d\n", iseg);
    ChebyModel_ReadFits(&cms->segments[iseg], f, &row0);
    //if ((ret=ChebyModel_ReadFits(&cms->segments[iseg], f, &row0)) != 0)
		//{
		//	printf ("ret %d\n", ret);
    //  return ret;
		//}
	}
  //printf("Got here f\n");
  return 0;
}

int ChebyModel_ReadFits(ChebyModel *cm, fitsfile *f, int *row0)
{
  int first = 1;
  char keyword[64], arg[64], junk[1024];
  int nx=-1, ny=-1, ix=0, iy;
  int ichar, nread;

	char **line;
	line = (char **)malloc(sizeof(char *));
	line[0] = (char *)malloc(sizeof(char)*1024);

	int colnum;
	long int nrows;
	int frow;
	int status = 0;
	int anynull = 0;
	char nval[]="NULL";

	if (fits_get_colnum(f,CASEINSEN,"PREDICT",&colnum,&status))
	{
		printf( "error while getting the column number\n" );
	}

	if ( fits_get_num_rows(f, &nrows, &status) )           // get the row number
	{
		printf( "error while getting the row number\n" );
	}

  cm->cheby.coeff=NULL;

	for (frow = (*row0); frow <= nrows; frow++)
	{
    //printf("Read %d\n", frow);
		fits_read_col(f, TSTRING, colnum, frow, 1, 1, nval, line, &anynull, &status);           // read the column
		//puts(line[0]);
    //printf("Read %s\n",line);
    if (sscanf(line[0], "%s", keyword)!=1)
      continue; // skip blank lines
    if (sscanf(line[0], "%s %s", keyword, arg)!=2)
      return -2;
    if (line[0][0]=='#')
      continue; // skip comment lines
    // check first line
    if (first && (strcasecmp(keyword, "ChebyModel")||strcasecmp(arg, "BEGIN")))
		{
      return -3;
		}
    // parse based on keyword
    if (!strcasecmp(keyword, "PSRNAME"))
      strcpy(cm->psrname, arg);
    else if (!strcasecmp(keyword, "SITENAME"))
      strcpy(cm->sitename, arg);
    else if (!strcasecmp(keyword, "TIME_RANGE"))
    {
      if (sscanf(line[0], "%*s %Lf %Lf", &cm->mjd_start, &cm->mjd_end)!=2)
				return -4;
    }
    else if (!strcasecmp(keyword, "FREQ_RANGE"))
    {
      if (sscanf(line[0], "%*s %Lf %Lf", &cm->freq_start, &cm->freq_end)!=2)
				return -5;
    }
    else if (!strcasecmp(keyword, "DISPERSION_CONSTANT"))
    {
      if (sscanf(arg, "%Lf", &cm->dispersion_constant)!=1)
				return -6;
    }
    else if (!strcasecmp(keyword, "NCOEFF_TIME"))
    {
      if (sscanf(arg, "%d", &nx)!=1) 
				return -7;
    }
    else if (!strcasecmp(keyword, "NCOEFF_FREQ"))
    {
      if (sscanf(arg, "%d", &ny)!=1)
				return -8;
    }
    else if (!strcasecmp(keyword, "COEFFS"))
    {
			//printf("Have nx = %d %d\n",nx,ny);
      if (cm->cheby.coeff==NULL) // first instance of COEFF keyword
			{
				if (nx < 0 && ny < 0) // oops, these should come first!
					return -8;
	
				ChebyModel_Init(cm, nx, ny);
				//printf("Initialised the memory\n");
      }
      if (ix >= nx)
				return -9; // too many coefficient lines!!

			sscanf(line[0], "%*s %n", &ichar);
      //printf("Line = %s %d %d\n",line,ny,ichar);
      if (ny<4) // All on one line
			{
				for (iy=0; iy < cm->cheby.ny; iy++)
				{
					if (sscanf(line[0]+ichar, "%Lf %n", &cm->cheby.coeff[iy*cm->cheby.nx+ix], &nread)!=1)
						return -10;
	      
					ichar += nread;
				}
			}
      else   // Code added by G. Hobbs for multiple lines in the predictor file
			{
				for (iy=0; iy < cm->cheby.ny; iy++)
				{
					if (sscanf(line[0]+ichar, "%Lf %n", &cm->cheby.coeff[iy*cm->cheby.nx+ix], &nread)!=1)
						return -10;
	      
					ichar += nread;
	      
					if ((iy+1)%3==0)
					{
						if (sscanf(line[0]+ichar, "%s", junk)==1)
							return -11; // excess stuff at end of line		  
		  
						ichar = 0;
						frow++;
						fits_read_col(f, TSTRING, colnum, frow, 1, 1, nval, line, &anynull, &status);           // read the column
						//puts(line[0]);
						//if (fgets(line, 1024, f)!=line)
						//	return -1;
					}
				}
			}
			if (sscanf(line[0]+ichar, "%s", junk)==1)
				return -11; // excess stuff at end of line
			ix++;
		}
    else if (!strcasecmp(keyword, "ChebyModel"))
    {
      if ((!first) && !strcasecmp(arg, "BEGIN"))
			{
				(*row0) = frow;
				return -12;
			}
      else if (!strcasecmp(arg, "END"))
      {
				if (cm->cheby.coeff==NULL || ix!=nx)
					return -13; // haven't read enough coefficients yet!!
				else
				{
					Cheby2D_Construct_x_Derivative(&cm->frequency_cheby, &cm->cheby);
				}
			}
    }
    else
      return -14; // unrecognized keyword!! 
    first = 0;    
  }

  return 0;
}

