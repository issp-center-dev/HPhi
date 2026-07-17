.. highlight:: none

保存量に関するパラメータ
--------------------------

*  ``nelec``

   **形式 :** 整数

   **説明 :** 伝導電子数を指定します。 ``model = "Fermion HubbardGC"``,
   ``"Spin"``, ``"SpinGC"`` のときには指定しないでください。

*  ``2Sz``

   **形式 :** 整数

   **説明 :** 全スピンのz 成分の2倍を指定します。
   ``model = "Fermion HubbardGC"``, ``"SpinGC"``
   のときには指定しないでください。

*  ``MomentumIndex``

   **形式 :** 整数

   **説明 :** 1次元並進対称性の運動量セクターを整数 :math:`m` で
   指定します。生成されるエキスパート入力には ``qptransidx.def`` と
   ``namelist.def`` の ``TransSym`` エントリが追加されます。長さ
   :math:`L` の chain で :math:`g` サイト並進したときの指標は
   :math:`\exp(-2\pi i m g/L)` です。 ``MomentumIndex`` は
   :math:`0 \le m < L` を満たす必要があります。

   現在は HPhi のスタンダードモードで ``lattice = "chain"``,
   ``phase0`` を使わない周期境界条件、かつ ``"Lanczos"`` または
   ``"CG"`` の計算手法にのみ対応します。対応する model は以下です。

   * ``model = "Spin"`` では ``2S = 1``、固定 ``2Sz``、かつ
     ``Jz = 0`` で磁場・一般相互作用・pair 項を含まない
     exchange-only のスピン相互作用に対応します。
   * ``model = "Fermion Hubbard"`` または ``"Hubbard"`` では、固定
     ``nelec`` と ``2Sz``、最近接ホッピング ``t``/``t0``、オンサイト
     ``U`` に対応します。化学ポテンシャルや磁場項
     (``mu``, ``h``, ``Gamma``, ``Gamma_y``)、遠距離ホッピング
     (``t'``, ``t''``)、非局所 Coulomb 項 (``V``/``CoulombInter``)、
     一般相互作用、pair 項は ``MomentumIndex`` と併用できません。

.. raw:: latex

   \newpage
