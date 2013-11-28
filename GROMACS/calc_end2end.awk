BEGIN {frame=0; ndna=0; nee=0; eedist=0; eeoc=0; 
FN="end2end_vs_t.dat"
printf("") > FN
} 
{ 
if ($0~".*t=.*") 
  {
    time=$NF;
  }
if ($3=="P" && $5==2+ndna*24) 
  {
    PA1[0]=$6; PA1[1]=$7; PA1[2]=$8;
  }; 
if ($3=="P" && $5==12+ndna*24) 
  {
    PB1[0]=$6; PB1[1]=$7; PB1[2]=$8;
  }; 
if ($3=="P" && $5==14+ndna*24) 
  {
    PB2[0]=$6; PB2[1]=$7; PB2[2]=$8;
  }; 
if ($3=="P" && $5==24+ndna*24) 
  {
    PA2[0]=$6; PA2[1]=$7; PA2[2]=$8; 
    for (i=0; i < 3; i++) CM1[i]=0.5*(PA1[i]+PA2[i]); 
    for (i=0; i < 3; i++) CM2[i]=0.5*(PB1[i]+PB2[i]); 
    eeone=sqrt((CM1[0]-CM2[0])^2+(CM1[1]-CM2[1])^2+(CM1[2]-CM2[2])^2);
    eeoc+=eeone;
    eedist+=eeone;
    ndna=ndna+1; 
    nee++; 
 }; 
  if ($1=="ENDMDL") 
    {
        printf("%.15G %.15G\n",time, eeoc/ndna/10) >> FN;
        frame++;
        eeoc=0;
	ndna=0;
    }
} 
END { print ("AVG END2END = ", eedist/nee/10, "( ", nee, ")")}
