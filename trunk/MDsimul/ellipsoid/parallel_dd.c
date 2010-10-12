#include<mdsimul.h>
#ifdef ED_PARALL_DD
/* ==== >>> routines to initialize structures <<< ==== */
unsigned long long int calc_cellnum(int ix, int iy, int iz)
{
  /* N.B. interleaving cells coordinates in binary represantion a unique cell  
     number can be calculated (see Miller and Luding J. Comp. Phys. 2003) */
  const unsigned int maxk=20;
  unsigned int ixl, iyl, izl, km;
  unsigned long long int im=0x1, bitx, bity, bitz, res;
  /* max number of cells is 2^(3*20)*/
  ixl = ix;
  iyl = iy;
  izl = iz;

  /* se n_bits(i) è il numero massimo di bits dell'unsigned int i
     allora di seguito si valuta max{n_bits(ix),n_bits(iy), n_bits(iz)}*/
  for (k=0; k < maxk; k++)
    {
      if (ixl==0x0 && iyl==0x0 && izl==0x0)
	{
	  km = k;
	  break;
	}
      ixl = ixl >> 1;
      iyl = iyl >> 1;
      izl = izl >> 1;
    }
  res = 0x0;
  for (k=0; k < km; k++)
    {
      bitx = ix & im;
      bity = iy & im;
      bitz = iz & im;
      res |= (bitz << (2*k+2)) | (bity << (2*k+1)) | (bitx << (2*k)); 
      im = im << 1;
    }
  return res;
}
int cellToRegion(unsigned long long int cell)
{
  /* implementare una ricerca dicotomica */
  int i;
  for (i=0; i < numOfProcs; i++)
    {
      if (cell < regionsArr[i])
	return i;
    }
}
void dd_init(void)
{
  /* queste dichiarazioni vanno poi rese globali */
  int *inRegion, dd_numreg, cc, iX, iY, iZ, ncp;
  /* 08/10/10: nel caso di codice parallelo è meglio non usare le multiple linked list poiche'
     il tutto si complica significativamente */
  /* numero totale di regioni (1 regione = 1 processo) */
  dd_numreg = num_of_processes;
  /* number of regions per process */
  ncp = (int) pow(2,(int)log2(cellsx*cellsy*cellsz/dd_numreg));
  if (OprogStatus.dd_numcells_per_proc == 0)
    {
      printf("ERROR: less than a cell per region, exiting...!\n");
      exit(-1);
    }
  /* inRegion[c] gives region associated to cell c */
  inRegion = malloc(sizeof(int)*cellsx*cellsy*cellsx);
  /* array contenente la fine di ogni regione nella sequenza associata alle celle
     tramite interleaving dei bit (vedi Sec. 3.1 S. Miller and S. Luding J. Comput. Phys. 193, 
     306-316 (2003) */
  regionsArr = malloc(sizeof(int)*dd_numreg);
  /* le regioni sono in generale dei parallelepipedi contenenti un numero intero 
     di celle ed individuati da 6 numeri (x1,x2), (y1,y2), (z1,z2) ossia x1i e x2i, all'inizio
     i parallelepipedi vengono scelti tutti uguali (a parte gli ultimi se il numero di regioni
     lungo un asse non è sottomultiplo del numero di celle. */
  for (i = 0; i < dd_numreg-1; i++)
    {
      /* regionsArr[i] contiene la fine della regione i-esima (in numero di celle), ossia
	 la prima cella che non gli appartiene. */
      regionsArr[i] = (i+1)*ncp;
    }
  /* all cells left assigned to last region */
  regionsArr[dd_numreg-1] = cellsx*cellsy*cellsz;
  if (regionsArr[dd_numreg-1] % 2 != 0)
    {
      printf("ERROR: total number of cell is %d, but it must be even in a parallel simulation\n", 
	     cellsx*cellsy*cellsz);
      exit(-1);
    }
  for (iZ = 0; iZ < cellsz; iZ++)
    for (iY = 0; iY < cellsy; iY++)
      for (iX = 0; iX < cellsx; iX++)
	{
	  cellnum = calc_cellnum(iX,iY,iZ);
	  inRegion[cellnum] = cellToRegion(cellnum);
	}
}
/* ==== >>> rollback <<< ==== */

#endif
