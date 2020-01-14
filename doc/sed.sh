#[s] for *.c
find ./src/*.c -type f -print0 | xargs -0 sed -i  -e 's/GC_child/GC/g'
#
find ./src/*.c -type f -print0 | xargs -0 sed -i  -e 's/child_omp/omp/g'
find ./src/*.c -type f -print0 | xargs -0 sed -i  -e 's/child_exchange/exchange/g'
find ./src/*.c -type f -print0 | xargs -0 sed -i  -e 's/child_pairhopp/pairhopp/g'
find ./src/*.c -type f -print0 | xargs -0 sed -i  -e 's/child_pairlift/pairlift/g'
find ./src/*.c -type f -print0 | xargs -0 sed -i  -e 's/child_CisAis/CisAis/g'
find ./src/*.c -type f -print0 | xargs -0 sed -i  -e 's/child_CisAit/CisAit/g'
find ./src/*.c -type f -print0 | xargs -0 sed -i  -e 's/child_CisAjt/CisAjt/g'
find ./src/*.c -type f -print0 | xargs -0 sed -i  -e 's/child_CisAisCisAis/CisAisCisAis/g'
find ./src/*.c -type f -print0 | xargs -0 sed -i  -e 's/child_CisAisCitAiu/CisAisCitAiu/g'
find ./src/*.c -type f -print0 | xargs -0 sed -i  -e 's/child_CisAitCiuAiv/CisAitCiuAiv/g'
find ./src/*.c -type f -print0 | xargs -0 sed -i  -e 's/child_CisAjtCkuAlv/CisAjtCkuAlv/g'
find ./src/*.c -type f -print0 | xargs -0 sed -i  -e 's/child_CisAjtCkuAku/CisAjtCkuAku/g'
find ./src/*.c -type f -print0 | xargs -0 sed -i  -e 's/child_general/general/g'
#
find ./src/*.c -type f -print0 | xargs -0 sed -i  -e 's/X_/child_/g'
#[e] for *.c

#[s] for *.h
find ./src/include/*.h -type f -print0 | xargs -0 sed -i  -e 's/GC_child/GC/g'
#
find ./src/include/*.h -type f -print0 | xargs -0 sed -i  -e 's/child_omp/omp/g'
find ./src/include/*.h -type f -print0 | xargs -0 sed -i  -e 's/child_exchange/exchange/g'
find ./src/include/*.h -type f -print0 | xargs -0 sed -i  -e 's/child_pairhopp/pairhopp/g'
find ./src/include/*.h -type f -print0 | xargs -0 sed -i  -e 's/child_pairlift/pairlift/g'
find ./src/include/*.h -type f -print0 | xargs -0 sed -i  -e 's/child_CisAis/CisAis/g'
find ./src/include/*.h -type f -print0 | xargs -0 sed -i  -e 's/child_CisAit/CisAit/g'
find ./src/include/*.h -type f -print0 | xargs -0 sed -i  -e 's/child_CisAjt/CisAjt/g'
find ./src/include/*.h -type f -print0 | xargs -0 sed -i  -e 's/child_CisAisCisAis/CisAisCisAis/g'
find ./src/include/*.h -type f -print0 | xargs -0 sed -i  -e 's/child_CisAisCitAiu/CisAisCitAiu/g'
find ./src/include/*.h -type f -print0 | xargs -0 sed -i  -e 's/child_CisAitCiuAiv/CisAitCiuAiv/g'
find ./src/include/*.h -type f -print0 | xargs -0 sed -i  -e 's/child_CisAjtCkuAlv/CisAjtCkuAlv/g'
find ./src/include/*.h -type f -print0 | xargs -0 sed -i  -e 's/child_CisAjtCkuAku/CisAjtCkuAku/g'
find ./src/include/*.h -type f -print0 | xargs -0 sed -i  -e 's/child_general/general/g'
#
find ./src/include/*.h -type f -print0 | xargs -0 sed -i  -e 's/X_/child_/g'
#[e] for *.h
