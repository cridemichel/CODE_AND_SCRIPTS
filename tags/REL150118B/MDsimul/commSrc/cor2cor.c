#include<mdsimul.h>
char inputFile[NAME_LENGTH];

/* ============================ >>> main <<< ============================= */
void main(int argc, char* argv[])
{
  int bakfd;

  /* input file is a the name of the backup file */
  if ( (bakfd = open(inputFile, O_RDONLY)) == -1)
    {
      printf("ERROR: Unable to open the coordinate file!\n");
      exit(-1);
    }

  readCoord();

  close(bakfd)

}
