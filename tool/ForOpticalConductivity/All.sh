#[s] execute HPhi
HPhi -s stan.in
#[e] execute HPhi

#[s] for negative omega
python3 OpticalSpectrum.py  1
HPhi -e Snamelist.def
cp -r output P
cp *.def ./P
#[e] for negative omega

#[s] for negative omega
python3 OpticalSpectrum.py  -1
HPhi -e Snamelist.def
cp -r output M
cp *.def ./M
#[e] for negative omega

#[s] sorting
sort -r -n P/zvo_DynamicalGreen.dat > tmp_P
cat        M/zvo_DynamicalGreen.dat > tmp_M
paste tmp_P tmp_M > optical.dat
#[e] sorting
