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
   * ``model = "SpinlessFermion"`` では固定 ``ncond`` と hopping-only
     の項に対応します。密度相互作用 (``V``/``CoulombInter``)、
     一般相互作用、pair 項はスタンダードモードの ``MomentumIndex``
     と併用できません。
   * ``model = "Fermion Hubbard"`` または ``"Hubbard"`` では、固定
     ``nelec`` と ``2Sz``、最近接ホッピング ``t``/``t0``、オンサイト
     ``U`` に対応します。化学ポテンシャルや磁場項
     (``mu``, ``h``, ``Gamma``, ``Gamma_y``)、遠距離ホッピング
     (``t'``, ``t''``)、非局所 Coulomb 項 (``V``/``CoulombInter``)、
     一般相互作用、pair 項は ``MomentumIndex`` と併用できません。

   各運動量点の最低エネルギーを求めるには、``MomentumIndex`` だけを
   変更し、:math:`m=0,\ldots,L-1` ごとにスタンダードモードの計算を
   個別に実行します。例えば6サイトのスピン鎖の :math:`m=1`
   セクターは次のように指定できます。

   .. code-block:: text

      L = 6
      model = Spin
      method = Lanczos
      lattice = chain
      Jx = 1.0
      Jy = 1.0
      Jz = 0.0
      2Sz = 0
      MomentumIndex = 1

   この入力を ``stan.in`` に保存した場合は ``HPhi -s stan.in`` を
   実行します。対応する運動量は :math:`k=2\pi m/L` で、最低エネルギーは
   ``output/zvo_energy.dat`` に出力されます。結果が上書きされないよう、
   ``MomentumIndex`` ごとに異なる作業ディレクトリで実行するか、次の計算前に
   出力を保存してください。

.. raw:: latex

   \newpage
