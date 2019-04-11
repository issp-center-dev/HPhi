動的グリーン関数の計算に関するパラメーター
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  ``CalcSpec``

   **形式 :** 文字列(\ ``"None"``, ``"Normal"``, ``"NoIteration"``,
   ``"Restart_out"``, ``"Restart_in"``,
   ``"Restart"``, ``"Scratch"``\ 。デフォルトは\ ``"None"``)

   **説明 :** 動的グリーン関数の計算に関する設定を行う。
   ``"None"``\ では動的グリーン関数を計算しない。
   ``"Normal"``\ では一から動的グリーン関数の計算を始める。
   ``"NoIteration"``\ では、前回の反復回数と同じところまで反復させる。
   このとき、ハミルトニアン-ベクトル積演算は行われないため、
   計算コストは非常に軽いが、十分な精度が出せない場合がある。
   ``"Restart_out"``\ では一から計算を始めて、
   反復が終了した時点で再計算用のデータをファイル出力する。
   ``"Restart_in"``\ では再計算用のデータをファイルから受け取り途中から計算を始める。
   ``"Restart"``\ では再計算用のデータをファイルから受け取り途中から計算を始め、
   反復が終了した時点で再計算用のデータをファイル出力する。
   ``"Scratch"``\ では、まず固有状態の計算が行われたのちに、
   続けて動的グリーン関数の計算が行われる。この場合は波動関数のファイル入出力は
   必要なくなる。
   スペクトル計算において使用される手法はパラメーター\ ``method``\ で指定されます。
   (``method="CG"``\ とした場合には
   付属している\ :math:`K\omega`\ ライブラリ [#]_ が呼び出され、
   シードスイッチ [#]_ 付きシフト双共役勾配法 [#]_ が使われます。

-  ``SpectrumType``

   **形式 :** 文字列(\ ``"SzSz"``, ``"S+S-"``, ``"Density"``, ``"up"``,
   ``"down"``, ``"SzSz_R"``, ``"S+S-_R"``, ``"Density_R"``, ``"up_R"``,
   ``"down_R"``\ のいずれか。デフォルトは\ ``"SzSz"``)

   **説明 :** 計算する動的グリーン関数の種類を指定する。
   逆格子空間での相関関数を直接計算するに場合には次の値を使用する。
   ``"SzSz"``\ では\ :math:`\langle {\hat S}_{z q} {\hat S}_{z q}\rangle`\ 、
   ``"S+S-"``\ では\ :math:`\langle {\hat S}^{+}_{q} {\hat S}^{-}_{q}\rangle`\ 、
   ``"Density"``\ では\ :math:`\langle {\hat n}_{q} {\hat n}_{q}\rangle`\ 、
   ``"up"``\ では\ :math:`\langle {\hat c}^{\dagger}_{q \uparrow} {\hat c}_{q \uparrow}\rangle`\ 
   ``"down"``\ では\ :math:`\langle {\hat c}^{\dagger}_{q \downarrow} {\hat c}_{q \downarrow}\rangle`
   。
   また、実空間での相関関数の計算では次の値を用いる。
   ``"SzSz_R"``\ では\ :math:`\langle {\hat S}_{z R} {\hat S}_{z 0}\rangle`\ 、
   ``"S+S-_R"``\ では\ :math:`\langle {\hat S}^{+}_{R} {\hat S}^{-}_{0}\rangle`\ 、
   ``"Density_R"``\ では\ :math:`\langle {\hat n}_{R} {\hat n}_{0}\rangle`\ 、
   ``"up_R"``\ では\ :math:`\langle {\hat c}^{\dagger}_{R \uparrow} {\hat c}_{0 \uparrow}\rangle`\ 、
   ``"down_R"``\ では\ :math:`\langle {\hat c}^{\dagger}_{R \downarrow} {\hat c}_{0 \downarrow}\rangle` 、
   ここで :math:`R` はすべてのサイト番号にわたる。
   この実空間での動的相関関数をフーリエ変換して逆格子空間全体での
   動的相関関数を得る方法については
   :ref:`相関関数のFourier変換ユーティリティー <fourier>` を参照。

-  ``SpectrumQW``, ``SpectrumQL``, ``SpectrumQH``

   **形式 :** 実数(デフォルトはともに\ ``0.0``)

   **説明 :** ``SpectrumType`` が ``"SzSz"``, ``"S+S-"``, ``"Density"``, ``"up"``,
   ``"down"`` のときのみ使用。
   計算する動的グリーン関数の波数を Fractional
   coordinateで指定する。 逆格子ベクトルは
   :numref:`fig_chap04_1_lattice`, :numref:`fig_chap04_1_honeycomb`,
   :numref:`fig_ladder`, :numref:`fig_kagome`
   に表されている格子ベクトルと対応するものとなる。

-  ``OmegaMin``

   **形式 :**
   実数(デフォルトは\ ``-LargeValue``\ :math:`\times`\ サイト数)

   **説明 :** 計算する動的グリーン関数の振動数の実部の下限。

-  ``OmegaMax``

   **形式 :**
   実数(デフォルトは\ ``LargeValue``\ :math:`\times`\ サイト数)

   **説明 :** 計算する動的グリーン関数の振動数の実部の上限。

-  ``OmegaIm``

   **形式 :** 実数(デフォルトは\ ``0.01*LargeValue``)

   **説明 :** 計算する動的グリーン関数の振動数の虚部。

-  ``NOmega``

   **形式 :** 正の整数(デフォルトは\ ``200``)

   **説明 :** 計算する動的グリーン関数の振動数のグリッド数。

.. [#] https://github.com/issp-center-dev/Komega.
.. [#] \S. Yamamoto, T. Sogabe, T. Hoshi, S.-L. Zhang, T. Fujiwara, Journal of the Physical Society of Japan **77**, 114713 (2008).
.. [#] \A. Frommer, Computing **70**, 87{109 (2003).

.. raw:: latex

   \newpage
