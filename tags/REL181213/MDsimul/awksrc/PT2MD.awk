#PROGRAMMA per riscalare misture monoatomiche e liquidi semplci (monoatomici)
# gawk -v tipo = {"pos","vel","for","blo"} -v blocco={1,2,3,...} -v S=<fattore>
# <file Cnf>
#NOTA: blocco ha senso solo se tipo = "blo"
BEGIN { prima = -1}

$1 == "parnum:" {  if (NF == 2) #se si tratta di liquido monoatomico ho un solo arg
			parnum = $2 
		else if (NF == 3)#se si tratta di una mistura ho A e B 
			parnum = $2 + $3}

$1 == "T:" { T = $2 }

$1 == "lambda0:" {  }#memorizza tutti i lambda0 in un vettore

$1 == "lambdat:" {  }#memorizza tutti i lambdat in un vettore

}
 
$1 == "@@@" { print $0}

$1 ~ ":" { print $0}

function scala()
{
	$1 = $1 * S	
	$2 = $2 * S
	$3 = $3 * S
	printf ("%.15G %.15G %.15G\n", $1,$2,$3)
}

function scalablocco()
{
if ((NR - prima) >= parnum*(blocco-1) && (NR - prima) < parnum*(blocco)){
		scala()
		}
	     else
		print $0
}
function scalapos() {
if ((NR - prima) >= 0 && (NR - prima) < parnum){
		scala()
		}
	     else
		print $0
}

function scalavel() {
if ((NR - prima) >= parnum && (NR - prima) < 2*parnum){
		scala()
		}
	     else
		print $0
}

function scalafor() {
if ((NR - prima) >= 2*parnum && (NR - prima) < 3*parnum){
		scala()
		}
	     else
		print $0
}

$1 !~ ":" && $1 != "@@@" { if (prima == -1) prima = NR; 
		if (tipo == "pos")
			scalapos()	    
		else if (tipo == "vel")
			scalavel()		 
		else if (tipo == "for")
			scalafor()
		else if (tipo == "blo")
			scalablocco()
		else
			print "ERRORE: TIPO SCONOSCIUTO!!!!"
	}
