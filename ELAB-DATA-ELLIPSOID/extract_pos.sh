tail -513 STORE/Store-$1-0 |head -256 |awk '{print $1,$2,$3}' > bbb.pos


paste aaa.pos bbb.pos |awk '{print $1,$4}' > rx.pos
paste aaa.pos bbb.pos |awk '{print $2,$5}' > ry.pos
paste aaa.pos bbb.pos |awk '{print $3,$6}' > rz.pos
