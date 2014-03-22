BEGIN {
  print "";
  #angstrom units (radius)
  Omass=16.0;
  Hmass=1.0;
  Cmass=12.0;
  Pmass=31.0;
  Nmass=14.0;
  msms=1;
  mass=0.0;
  comx= 0;
  comy= 0;
  comz= 0;
  hydration=0;
}
{
 format="%7.5G    %7.5G    %7.5G  %G";
 got=0;

 rx = $7;
 ry = $8;
 rz = $9;
 
 #rx = rx - comx;
 #ry = ry - comy;
 #rz = rz - comz; 
 if ($1=="HETATM")
  {
    if (hydration==1)
     {
       mass+=18; comx+=18*rx; comy+=18*ry; comz+=18*rz;
       printf(format,rx,ry,rz,18);
       if (msms==0)
        print $3;
       else 
         print "";
      }
 }
 else 
 {
 if ($3~"^H.*") {mass+=Hmass; comx+=Hmass*rx; comy+=Hmass*ry; comz+=Hmass*rz;  printf(format, rx,ry,rz,Hmass); got=1};
 if ($3~"^O.*") {mass+= Omass; comx+=Omass*rx; comy+=Omass*ry; comz+=Omass*rz;  printf(format, rx,ry,rz,Omass)  ; got=1};
 if ($3~"^P.*") {mass+=Pmass; comx+=Pmass*rx; comy+=Pmass*ry; comz+=Pmass*rz;  printf(format, rx,ry,rz,Pmass)  ; got=1};
 if ($3~"^C.*") {mass+=Cmass; comx+=Cmass*rx; comy+=Cmass*ry; comz+=Cmass*rz;  printf(format, rx,ry,rz,Cmass)  ; got=1};
 if ($3~"^N.*") {mass+=Nmass; comx+=Nmass*rx; comy+=Nmass*ry; comz+=Nmass*rz;  printf(format, rx,ry,rz,Nmass)  ; got=1};
 if ($3~"^1H.*"){mass+= Hmass; comx+=Hmass*rx; comy+=Hmass*ry; comz+=Hmass*rz;  printf(format, rx,ry,rz,Hmass)  ; got=1};
 if ($3~"^2H.*"){mass+= Hmass; comx+=Hmass*rx; comy+=Hmass*ry; comz+=Hmass*rz;  printf(format, rx,ry,rz,Hmass)  ; got=1};
 if (got==1)
  {
    if (msms==0)
      print $3;
    else 
      print "";
  }
 }
}

END {
printf("CoM - Mtot = %f %f %f %f\n", comx/mass, comy/mass, comz/mass, mass);
}
