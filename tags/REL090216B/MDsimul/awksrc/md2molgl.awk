# si puo' usare per sistemi misture biamtomiche e liquidi monoatomici 
BEGIN {livgrA = 50; livgrB = 100;
       print ".Vol: 256\n.numAt: 1"; prima = -1}

function printpos() {
if (MIX)
{
if ((NR - prima) >= 0 && (NR - prima) < parnumA){  
   
   print ($0 " @" raggioA " $" livgrA)
	}
else if ((NR - prima) >= parnumA && (NR - prima) < parnumB) {
   print ($0 " @" raggioB " $" livgrB)

   }
}
else {
if ((NR - prima) >= 0 && (NR - prima) < parnumA){  
   print ($0 " @" raggio)
}
}

$1 ~ "sigma" { if (NF == 2)
                {  raggio = $2
                }
               else
                 { raggioA = $2; raggioB = $5}
                }
$1 ~ "parnum" { if (NF == 2) #se si tratta di liquido monoatomico ho un solo arg
			{
                          MIX = 0;
                          parnum = $2
			    
                         } 
		else if (NF == 3)#se si tratta di una mistura ho A e B 
			{
                            MIX = 1;
                            parnumA = $2; 
			    parnumB = $3;
			    parnum  = parnumA + parnumB;
			}

$1 !~ "*:" && $1 != "@@@" { if (prima == -1) prima = NR;
                             printpos() }
