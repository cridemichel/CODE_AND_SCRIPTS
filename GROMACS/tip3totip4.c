#include <stdlib.h>
#include <stdio.h>
char resname[256], atomtype[256], line[4096], nresstr[256], atomnumbstr[256]; 
int nres, atomnumb, atno=0, i;
const char mwtype[]="MW";
double rx, ry, rz, vx, vy, vz, Lx, Ly, Lz;
int main(int argc, char** argv)
{
 FILE *fi, *fo;
 fi=fopen(argv[1], "r");
 fo=fopen(argv[2], "w+");
 
 fscanf(fi, "%[^\n]\n", line);
 fprintf(fo, "%s\n", line);
 fscanf(fi, "%d%c", &atomnumb, line);
 fprintf(fo, "%d\n", atomnumb);
 
 //printf("line=%s line[0]=%c atomtype=<%s>\n", line, line[0], atomtype);
 while (!feof(fi))
 { 
       fscanf(fi, "%69c", line);
  
 
	 //printf("line=*%s* line[0]=%c atomtype=<%s>\n", line, line[0], atomtype);
 if (sscanf(line, "%5c%5c%5c%5c%lf%lf%lf%lf%lf%lf", nresstr, resname, 
	atomtype, atomnumbstr, &rx, &ry, &rz, &vx, &vy, &vz) < 10)
 {	
    sscanf(line, "%lf %lf %lf\n", &Lx, &Ly, &Lz);
    fprintf(fo, "%f %f %f\n", Lx, Ly, Lz);
 }	
 else {
	 nres = atoi(nresstr);
         atomnumb=atoi(atomnumbstr);
	 atno++; 
	 fprintf(fo, "%5i%5s%5s%5i%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n", nres, resname, atomtype, atno, rx, ry, rz, vx, vy, vz);
	 if (!strcmp(atomtype, "  HW2"))
	 {
		 atno++;
		 fprintf(fo, "%5i%5s%5s%5i%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n", nres, resname, mwtype, atno, rx, ry, rz, vx, vy, vz);
	 }
 }
      for (i=0; i < 69; i++) line[i]='\n';
 }
 fclose(fi);
 fclose(fo);
 return 0; 
} 
