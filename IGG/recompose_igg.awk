{
   if (NF==13 && $13==0)
    {
      xref=$1;
      yref=$2;
      print $0;
    }
   else if (NF==13 && $13 > 0 && $13 < 4)
    {	delx=$1-xref;
        dely=$2-yref;
        if (delx < 0 && -delx > Lbox*0.5)
           printf("%.15G ", $1+Lbox);	
        else if (delx > Lbox*0.5)
           printf("%.15G ", $1-Lbox);	
        else 
           printf("%.15G ", $1);

        if (dely < 0 && -dely > Lbox*0.5)
           printf("%.15G ", $2+Lbox);	
        else if (dely > Lbox*0.5)
           printf("%.15G ", $2-Lbox);	
        else
           printf("%.15G ", $2);
        print($3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13);
    }
   else 
     print $0;
}
