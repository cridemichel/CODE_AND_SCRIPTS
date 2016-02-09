#define MAXINPS 255
#define MAXDAT 1000000
char outFile[512];
char inputFiles[MAXINPS][512];
int numinf=0, numx = 0;
/* ====================== >>> invalidArg <<< ===============================*/
void invalidArg(void)
{
  printf("Usage: mean [-o <mean_file>] <file1> <file2> ...\n");
  exit(-1);
}

/* =============================== >>> args <<< ============================ */
void args(int argc, char* argv[])
{
  int i, j;
  j = 0;
  for(i=1; i<argc; ++i)
    {
      if (!strcmp(argv[i],"-o"))
	{
	  if (i+1 == argc)
	    invalildArg();
	  ++i;
	  strcpy(outFile, argv[i]);
	  continue;
	}
      strcpy(inputFiles[j], argv[i]);
      ++j;
      if (j > MAXINPS) 
	{
	  printf("Sorry, maximum nnumber of files to average exceede!\n");
	  exit(-1);
	}
    }
  numinf = j;
}

/* =============================== >>> main <<< =============================*/
void main(int argc, char* argv[])
{
  double* inpdatx, valy;
  int* inpdaty;
  int i, valx;
  char line[1024];
  FILE* offs;
  FILE* iffs[MAXINPS];
  args(argc, argv);
 
  inpdatx = malloc(numinps * sizeof(double));
  inpdaty = malloc(numinps * sizeof(int)); 
  for(i = 0; i < MAXDAT; ++i)
    {
      inpdaty = 0.0;
    }

  offs = fopen(outFile, "w");
  for(i=0; i < numinf; ++i)
    {
      iffsp[i] = fopen(inputFiles, "r");
    }
  
  for(i=0; i<numinf; ++i)
    {
      j = 0;
      while(!feof(iffs[i]))
	{
	  fscanf(iffs[i],"%[^\n]\n", line);
	  if (!strcmp(line, "&"))/* skip if line = "&" */ 
	    continue;
	  sscanf(line, "%f %f\n", valx, valy);
	  if (inpdatx[j] == 0) inpdatx[j] = valx;
	  else
	    if (inpdatx[j] != valx)
	      {
		printf("Incompatible data files!\n");
		exit(-1);
	      }
	  inpdaty[j] += valy;
	  ++j;
	}
      numx = j;
    }

  for (i=0; i < numx; ++i)
    inpdaty[i] /= numinps;

  /* Write mean to disk */
  for (i=0; i< numx; ++i)
    {
      fprint(offs,"%.6f %.6f\n", inpdatx[i], inpdaty[i]);
    }

  printf("File with mean wroten\n");

  /* Close all streams */
  foclose(offs);
  for(i=0; i<numinf; ++i)
    {
      fclose(iffs[i]);
    }
  
}
