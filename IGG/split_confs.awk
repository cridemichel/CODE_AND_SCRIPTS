BEGIN { numant=0; iggnum=-1; NRB=-1; NRE=-1; at=0;}
{ 
if (at < 2)
{ 
  if ($1=="250" && $2=="250" && $3=="250" && $4=="250" && $5=="2050")
    print ("1 1 1 1 _FIXNUMANT_");
  else if ($1=="parnum:") 
    print ("parnum: _FIXPARNUM_");
  else print ($0);
} 
if (NR >= NRB && NR < NRB+Nigg && NF==13)
  linea[$13]=$0;
if (iggnum==-1 && NR >= NRB && NR < NRB+Nigg && $1 > xm && $1 < xM && $2 > ym && $2 < yM && $3+Lbox*0.5 < dl) 
{
  if ($13!=0)
   {
     for (kk=0; kk < $13; kk++)
       print(linea[kk]);
  }	
  iggnum=NR-$13;
  print ("TYPE=", $13) > "/dev/stderr" 
}
#print out a single igg coordinates
if (iggnum !=-1 && NR >= iggnum && NR < iggnum+4)
  {
     xc = $1-xm-dl*0.5;
     yc = $2-ym-dl*0.5;
     zc = $3;#+Lbox*0.5;
     # apply PBC 
     if (xc > dl*0.5) xc = xc - dl;
     if (xc < -dl*0.5) xc = xc + dl;
     if (yc > dl*0.5) yc = yc - dl;
     if (yc < -dl*0.5) yc = yc + dl;
     print (xc,yc,zc,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13); 
  }
#print out igg velocities
if (iggnum!=-1 && NR >= iggnum+np && NR < iggnum+np+4)
  print $0;
#print out antigens
if (NR-NRB >= Nigg && NR < NRE && $1 > xm && $1 < xM && $2 > ym && $2 < yM)
 {
   numant++;
   print $0;
 }
##if (NR-NRB >= Nigg+np && NR < NRE+np)
##  print $0;
if ($0=="@@@") 
   at++; 
if (at==2 && NRB==-1) 
   {NRB=NR+1; NRE=NR+np+1;};  
}
END {
   print ("NRB= ", NRB, "NRE=", NRE, "iggnum=", iggnum) > "/dev/stderr"; 
   print("xm=", xm, " xM=", xM, " ym=",ym, " yM=",yM)> "/dev/stderr"; 
    #print out box size
   for (kk=0; kk < numant; kk++)
     print ("0 0 0 0 0 0");
   printf ("%d\n",numant) > "./_NUMANT_";
   print(dl,dl,dl);
   if (iggnum==-1)
    {print "FAILED"};
 }
