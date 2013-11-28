BEGIN {frame=0; ndna=0; cc=0; diam=0; npairs=0; 
Pradius=0.18;
FN="end2end_vs_t.dat"
printf("") > FN
} 
{ 
if ($0~".*t=.*") 
  {
    time=$NF;
  }
  if ($3=="P") 
  { 
    if ($5-ndna*24 <= 12)
     {
	PA[$5-2,0]=$6; PA[$5-2,1]=$7; PA[$5-2,2]=$8;
        cc++;
	npairs++;
     }
   else
     {
	PB[24-$5,0]=$6; PB[24-$5,1]=$7; PB[24-$5,2]=$8;
     }
  if ($5==24+ndna*24) 
    {
      #printf("npairs=%d 5=%d  %d\n", npairs, $5, 24+ndna*24);
      for (i=0; i < npairs; i++) 
         diam+=sqrt((PA[i,0]-PB[i,0])^2+(PA[i,1]-PB[i,1])^2+(PA[i,2]-PB[i,2])^2);
      ndna=ndna+1; 
      npairs=0; 
     }; 
  }; 
 
  if ($1=="ENDMDL") 
    {
        frame++;
	ndna=0;
    }
} 
END { print ("AVG DIAM = ", diam/cc/10, "( ", cc, ")")}
