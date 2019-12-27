動的グリーン関数の計算に関するパラメータ
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  ``CalcSpec``

   **形式 :** 文字列(\ ``"None"``, ``"Normal"``, ``"NoIteration"``,
   ``"Restart_out"``, ``"Restart_in"``,
   ``"Restart"``\ 。デフォルトは\ ``"None"``)

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
   スペクトル計算において使用される手法はパラメーター\ ``method``\ で指定されます。
   (``method="CG"``\ とした場合には
   付属している\ :math:`K\omega`\ ライブラリ [#]_ が呼び出され、
   シードスイッチ [#]_ 付きシフト双共役勾配法 [#]_ が使われます。

-  ``SpectrumType``

   **形式 :** 文字列(\ ``"SzSz"``, ``"S+S-"``, ``"Density"``, ``"up"``,
   ``"down"``\ のいずれか。デフォルトは\ ``"SzSz"``)

   **説明 :** 計算する動的グリーン関数の種類を指定する。
   ``"SzSz"``\ では\ :math:`\langle {S}^z_{- \bf q} {S}^z_{\bf q}\rangle`\ 、
   ``"S+S-"``\ では\ :math:`\langle S^{+}_{- \bf q} {S}^{-}_{\bf q}\rangle`\ 、
   ``"Density"``\ では\ :math:`\langle {n}_{- \bf q} {n}_{\bf q}\rangle`\ 、
   ``"up"``\ では\ :math:`\langle {c}^{\dagger}_{{\bf q} \uparrow} {c}_{{\bf q} \uparrow}\rangle`\ 、
   ``"down"``\ では\ :math:`\langle {c}^{\dagger}_{{\bf q} \downarrow} {c}_{{\bf q} \downarrow}\rangle`
   となる。

-  ``SpectrumQW``, ``SpectrumQL``

   **形式 :** 実数(デフォルトはともに\ ``0.0``)

   **説明 :** 計算する動的グリーン関数の波数を Fractional
   coordinateで指定する。 逆格子ベクトルは
   :numref:`fig_chap04_1_lattice`, :numref:`fig_chap04_1_honeycomb`,
   :numref:`fig_kagome`, :numref:`fig_ladder`
   に表されている格子ベクトルと対応するものとなる。

-  ``OmegaOrg``

   **形式 :** 実数(デフォルトは\ ``0.0``)

   **説明 :**
   動的グリーン関数を計算する際の振動数\ :math:`\omega`\ の実部の原点。

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
.. [#] \A. Frommer, Computing **70**, 87-109 (2003).

.. raw:: latex

   \newpage
