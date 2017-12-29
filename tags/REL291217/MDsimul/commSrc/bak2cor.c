#include<mdsimul.h>
char inputFile[NAME_LENGTH];
char outFile[NAME_LENGTH];

/* ============================ >>> help <<< ============================= */
void help_b2c(void)
{
  printf("Syntax: bak2cor -o <output file> <coordinate file>\n");
}

/* ========================= >>> invalArg <<< ============================== */
void invalArg_b2c(void)
{ 
  /* DESCRIPTION:
     call this function when during parsing something goes wrong */
  help_b2c();
  printf("Invalid argument!\n");
  exit(-1);
}

/* =============================== >>> args <<< =============================*/
void args_b2c(int argc,char **argv)
{
  int i; /* fictitiuos counter */
  
  /* ==== set defaults =====*/
  
  strcpy(outFile, "scaled.cor"); 
  /* this is the defaukt output filename */

  /* =======================*/
  i = 1; 
  printf("argc: %d i: %d\n", argc, i);
  if (argc == 1) 
    invalArg_b2c();

  while ( i < (argc - 2) ) 
    { 
      if (!strcmp(argv[i], "-h"))
	{
	  help_b2c();
	  exit(-1);
	}
      else if (!strcmp(argv[i], "-o")) 
	{
	  if ( i+1 >= (argc - 1) )    /* argv[argc - 1] is the last arg, 
					 that mnst be the input file !!!! */ 
	    {
	      invalArg_b2c(); 
	    }
	  ++i;
	  strcpy(outFile, argv[i]);
	} 
      else                            /* option not existant */ 
	{
	  invalArg_b2c();
	}
    }
  /* and now take the input file */
  
  /* there is another arg ?  if N => ERROR */
  if (++i == argc)  
    {
      invalArg_b2c();
    }
  /* last arg that is input file name */
  strcpy(inputFile, argv[i]);
}

/*  ======================= >>> Function prototypes <<< =====================*/
void args(int argc, char **argv);
void loadCor(int cfd);
void saveCor(char* fileName);

/* ============================ >>> main <<< ============================= */
void main(int argc, char* argv[])
{
  int bakfd, corfd, bytes;
  char *buf;

  args_b2c(argc, argv);

  printf("Input File: %s\n", inputFile);
  /* input file is a the name of the backup file */
  if ( (bakfd = open(inputFile, O_RDONLY)) == -1)
    {
      printf("ERROR: Unable to open the backup file!\n");
      exit(-1);
    }

  if ( (corfd = open(outFile, O_CREAT | O_WRONLY, 0666)) == -1)
    {
      printf("ERROR: Unable to open the coordinate file!\n");
      exit(-1);
    }

  
  printf("Loading coordinates...\n");

  buf = malloc(sizeof(struct params));
  /* Read the Oparams structure */
  if (read(bakfd, buf, sizeof(struct params)) == -1)
    {
      perror("Reading Oparams header");
      exit(-1);
    }
  
  write(corfd, buf, sizeof(struct params));
  
  /* Skip the OprogStatus structure */
  lseek(bakfd, sizeof(struct progStatus), SEEK_CUR);
  
  free(buf);
  buf = malloc(1024);
  do
    {
      /* Read */
      if ( (bytes = read(bakfd, buf, 1024)) == -1)
	{
	  perror("Reading coordinates");
	  exit(-1);
	}
      /* Write */
      if (write(corfd, buf, bytes) == -1)
	{
	  perror("Writing coordinates");
	  exit(-1);
	}
    }
  while(bytes > 0); /* While not end of file */
  close(bakfd);
  close(corfd);
  printf("Done\n");
}
