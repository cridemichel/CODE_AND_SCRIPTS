#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#define PEP 300
#define nmaxfiles 16
char nomifiles[nmaxfiles][132];
double temp[nmaxfiles], PE[nmaxfiles][PEP], EE[PEP];
double PEij[nmaxfiles][PEP];
double bot[5];
const double bc1 = 14.0/45.0, bc2 = 64.0/45.0, bc3 = 24.0/45.0;
const Nm = 256;

/* =========================== >>> BodeTerm <<< ============================*/
double BodeTermB(double* fi)
{
  return (bc1 * fi[0] + bc2 * fi[1] + bc3 * fi[2] + bc2 * fi[3] +
	  bc1 * fi[4]);
}

 
int main(int argc, char** argv)
{
  FILE *inpf, *of;
  int ss, maxiE, iE, i, pe, ind, aggiungi, nn, in, nfiles, maxt;
  double norm, dblE, betai, betaj, maxE, maxPE, dE, ee, beta0;

  inpf = fopen(argv[1],"r");
  if (inpf == NULL)
    {
      printf("errore nell'apertura del file %s\n", argv[1]);
      exit(-1);
    }

  ind = 0;
  /* NOTA: il file con la lista di files deve essere del tipo:
     <nome_file_PE> <temperatura>  
   */
  while (1)
    {
      if (fscanf(inpf, "%s %lf\n", nomifiles[ind], &temp[ind]) < 2)
	break;
      ind++;
    }
  fclose(inpf); 
  nfiles = ind;
  /* legge tutte le PE da riscalare sulla prima */
  for (ind = 0; ind < nfiles; ind++)
    {
      inpf = fopen(nomifiles[ind], "r");
      if (inpf == NULL)
	{
	  printf("[%d] errore nell'apertura del file %s\n", ind, nomifiles[ind]);
	  exit(-1);
	}
      fprintf(stderr, "Leggo dal file %s\n", nomifiles[ind]);
      for (i = 0; i < PEP; i++)
	{
	  if (fscanf(inpf, "%lf %d", &ee, &pe) < 2)
	    break;
	  //printf("(%f,%d)\n", ee, pe);
	  PE[ind][i] = (double) pe;
	  EE[i] = ee;
	}
      fclose(inpf);
    }
 
  /* Questa e' il livello al di sotto del quale i valori delle
     distribuzioni di energia vengono tagliati */
  /* Se non c'e' abbastanza statistica esci */
  maxPE = 0.0;
  maxiE = 0;
  dE = fabs(EE[1] - EE[0]);

  /* Determina il massimo della curva di riferimento e usa il valore 
     dell'energia del massimo come riferimento per l'energia */
  for (iE = 0; iE < PEP; iE++)
    {
      if (PE[0][iE] > maxPE)
	{
	  maxiE = iE;
	  maxPE = PE[0][iE];
	}
    }

  maxE = EE[maxiE];
  betaj = 1.0/temp[0];
  for (ss = 0; ss < nfiles; ss++)
    {
      /* Usa per il check dell'equilibratura i quattro 
	 sistemi a piu' bassa temperatura */
      betai = 1.0 / temp[ss];
      /* Fattore di normalizzazione */
      for (iE = 0; iE < PEP; iE++)
	{
	  /* Calcola valori non senza normalizzazione */
	  dblE = ((double)Nm) * EE[iE];
	  dblE -= ((double)Nm) * maxE;
	  PEij[ss][iE] = ((double) PE[ss][iE])*
	    exp((betai-betaj)*dblE);
	}
      /* Calcola il fattore di normalizzazione integrando con il metodo
	 di Bode (ved. Numerical Recipe) */
      norm = 0.0;
      for (iE = 0; iE < PEP - 4; iE = iE + 4)
	{
	  for (i = 0; i < 5; i++)
	    {
	      /* Notare che l'energia e' opportunamente shiftata per 
		 avere numeri piu' piccoli a causa dell'esponeziale 
		 (si rischia l'overflow!!!) */
	      bot[i] = PEij[ss][iE+i];
	    }
	  /*
	    dblE = ((double)Nm) * (EN_MIN + ((double)iE+1 - maxiE) 
	    * (EN_MAX - EN_MIN) / 
	    POINTS));
	    PEij[ss][iE+1] = ((double) OprogStatus.PE[ss][iE+1])*
	    exp(beta0*(lambdai-lambdaj)*dblE);
	    norm += PEij[ss][iE] + PEij[ss][iE+1]; */
	  norm += BodeTermB(bot);
	}
      /* moltiplica per il deltaE */
      norm *= dE;
      
      for (iE = 0; iE < PEP - 1; iE++)
	{
	  if (norm != 0)
	    {
	      /* Normalizza! */
	      PEij[ss][iE] /= norm;
	    }
	}
    }
  /* stampa su stdout le nfiles distribuzioni riscalate */
  for (ind = 0; ind < nfiles; ind++)
    {
      for (i = 0; i < PEP; i++)
	{
	  fprintf(stdout, "%.8f %.15G\n", EE[i], PEij[ind][i]);
	}
      fprintf(stdout, "&\n");
    }
}

