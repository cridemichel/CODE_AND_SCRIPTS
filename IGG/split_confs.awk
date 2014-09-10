BEGIN { iggnum=0; }
{ 
if (NR > NRB && NR <= NRE && $1 > xm && $1 < xM && $2 > ym && $2 < yM) 
  iggnum=NR;
if (NR==iggnum+NRB)
  print ($1-xm-dl*0.5,$2-ym-dl*0.5,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13); 
if (NR==iggnum+NP+NRB)
  print $0;
if (NR >= NP+NRB && NR <=NRB+2*NP)
  print $0
if ($0=="@@@") 
   at++; 
if (at==2) 
   {NRB=NR; NRE=NR+NP};  
}
END
{
   print(dl,dl,dl);
}
