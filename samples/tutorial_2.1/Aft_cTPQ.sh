num_sample=10
for bs in  5
do
  python3 BS_TPQ.py  $num_sample $bs 
  python3 Ext.py  BS_MaxBS$bs.dat
done
