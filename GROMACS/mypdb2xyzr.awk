BEGIN {
  #angstrom units (radius)
  Orad=1.52;
  Hrad=1.2;
  Crad=1.7;
  Prad=1.8;
  Nrad=1.55;
 msms=1;
}
{
 got=0;
 if ($3~"^H.*") {printf("  %s   %s    %s %G ",$6,$7,$8,Hrad); got=1};
 if ($3~"^O.*") {printf("  %s   %s    %s %G ",$6,$7,$8,Orad)  ; got=1};
 if ($3~"^P.*") {printf("  %s   %s    %s %G ",$6,$7,$8,Prad)  ; got=1};
 if ($3~"^C.*") {printf("  %s   %s    %s %G ",$6,$7,$8,Crad)  ; got=1};
 if ($3~"^N.*") {printf("  %s   %s    %s %G ",$6,$7,$8,Nrad)  ; got=1};
 if (got==1)
  {
    if (msms==0)
      print $3;
    else 
      print "";
  }
}
