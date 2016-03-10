#include<stdlib.h>
#include"mdsimul.h"
#include<string.h>
void main(int argc, char** argv)
{
  char command[NAME_LENGTH];

  if (argc == 1) 
    {
      printf("Too few arguments\n");
      printf("usage: log md/chk to see mdsimul log file or\n");
      printf("check_status log file respectively\n");
      exit(-1);
    }
  if (!strcmp(argv[1],"md"))
      {
	strcpy(command, "less ");
	strcat(command, MD_HD_TMP);
	strcat(command, MDS_LOG);
	system(command);
      }
   else if (!strcmp("chk",argv[1]))
     {
	strcpy(command, "less ");
	strcat(command, MD_HD_TMP);
	strcat(command, CHK_LOG);
	system(command);
     }
  else 
    {
      printf("Invalig argument!\n");
    }
}
