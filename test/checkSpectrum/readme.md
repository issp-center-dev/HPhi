Lanczosスペクトル計算  
1. Lanczosの固有ベクトルをはき出す。  
    * calcmod.defで`OutputEigenVec 1`に設定する。  
2.  励起演算子用入力ファイルを作成し、namelist.defで指定する。  
3. 入力ベクトルを namelist.defで指定する(\_rankxxの前までを指定)。    
    * SpectrumVec zvo_eigenvec_0  
4. 励起計算モードに設定する。  
    * calcmod.defでCalcSpecの設定をする。  
    * `CalcSpec 1`に設定する(CalcSpecのその他オプションについてはマニュアルのキーワードの説明を参考にしてください)。  

TPQスペクトル計算  
1. TPQの再計算用ベクトルをはき出す。
    * calcmod.defで ` Restart 1 ` に設定する。  
2. 励起演算子用入力ファイルを作成し、namelist.defで指定する。  
3. 入力ベクトルを namelist.defで指定する。TPQの場合はsetの箇所まで指定する。  
    * SpectrumVec tmp_vec_set0  
4. 励起計算モードに設定する。  
    * calcmod.defでCalcSpecの設定をする。  
    * `CalcSpec 1`に設定する(CalcSpecのその他オプションについてはマニュアルのキーワードの説明を参考にしてください)。

全対角化スペクトル計算

 以下全てシングルプロセスで行う

1. Lanczosの固有ベクトルをはき出す。  
    * calcmod.defで`OutputEigenVec 1`に設定する。  
2.  励起演算子用入力ファイルを作成し、namelist.defで指定する。  
3. 入力ベクトルを namelist.defで指定する(\_rankxxの前までを指定)。
    * SpectrumVec zvo_eigenvec_0  
4. 励起計算モードに設定する。  
    * calcmod.defでCalcSpecの設定をする。  
    * `CalcSpec 1`に設定する(CalcSpecのその他オプションについてはマニュアルのキーワードの説明を参考にしてください)。
5. calcmodファイルでcalcTypeを 2 (全対角化)とする。