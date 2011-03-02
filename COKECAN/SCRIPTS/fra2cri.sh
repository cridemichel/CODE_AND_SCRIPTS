cat $1 | awk '{if ($17=="C[red]") print ($1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,"0")}' > _pos_
cat $1 | awk '{if ($17=="C[red]") print ("0 0 0 0 0 0")}' > _vel_
VOL="30"
cat header.cnf _pos_ _vel_ $VOL
rm -f _pos _vel_

