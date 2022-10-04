#[s] execute HPhi
$1 -s stan1.in
#[e] execute HPhi

#[s] for negative omega
python3 OpticalSpectrum.py  1
$1 -e Snamelist.def
cp -r output DG_P
cp *.def ./DG_P
#[e] for negative omega

#[s] for negative omega
python3 OpticalSpectrum.py  -1
$1 -e Snamelist.def
cp -r output DG_M
cp *.def ./DG_M
#[e] for negative omega

#[s] sorting
sort -r -n DG_P/zvo_DynamicalGreen.dat > tmp_P
cat        DG_M/zvo_DynamicalGreen.dat > tmp_M
paste tmp_P tmp_M > optical.dat
#[e] sorting
