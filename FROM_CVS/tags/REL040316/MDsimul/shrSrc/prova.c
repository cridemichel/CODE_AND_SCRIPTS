#include<share.h>
#include<stdio.h>
#include<stdlib.h>
void main(void)
{ char *prstr,*prstrB,*prstrC;
prstr = getshmem(10000);
prstrB = getshmem(10000);
prstrC = getshmem(10000);
prstr[0]='a';
printf("assegnato il primo byte di memoria condivisa:%c\n",prstr[0]);
freeshmem(prstr);
freeshmem(prstrB);
freeshmem(prstrC);
}



