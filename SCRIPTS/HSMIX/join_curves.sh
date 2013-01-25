i=0
while [ $[i] -lt 100 ]
do
cat Fqs-$i | head -n 20 > aaa
cat aaa Fqs-lin-$i > Fqs-tot-$i
ii=`echo $i | awk '{printf("%.3d",$0)}'` 
cat N-sqt.k=$ii | head -n 20 > aaa
cat aaa N-sqt-lin.k=$ii > N-sqt-tot.k=$ii
rm aaa
i=$[$i+1]
done
