.. highlight:: none

計算の種類に関する必須パラメーター
----------------------------------

*  ``model``

   **形式 :** 文字列(\ ``"Fermion Hubbard"``, ``"Spin"``,
   ``"Kondo Lattice"``, ``"Fermion HubbardGC"``, ``"SpinGC"``,
   ``"Kondo LatticeGC"``, ``"SpinGCCMA"``\ のいずれか) [#]_

   **説明 :** 計算対象の模型を指定します。上記の文字列はそれぞれ
   カノニカル集団のフェルミ粒子Hubbard模型

   .. math::
      :label: fml4_1_hubbard

      \mathcal H = -\mu \sum_{i \sigma} c^\dagger_{i \sigma} c_{i \sigma} 
      - \sum_{i \neq j \sigma} t_{i j} c^\dagger_{i \sigma} c_{j \sigma} 
      + \sum_{i} U n_{i \uparrow} n_{i \downarrow}
      + \sum_{i \neq j} V_{i j} n_{i} n_{j},

   同じくカノニカル集団のスピン模型
   (:math:`\{\sigma_1, \sigma_2\}={x, y, z}`)

   .. math::
      :label: fml4_1_spin

      \mathcal H &= -h \sum_{i} S_{i z} - \Gamma \sum_{i} S_{i x} + D \sum_{i} S_{i z} S_{i z}
      \nonumber \\
      &+ \sum_{i j, \sigma_1}J_{i j \sigma_1} S_{i \sigma_1} S_{j \sigma_1}+ \sum_{i j, \sigma_1 \neq \sigma_2} J_{i j \sigma_1 \sigma_2} S_{i \sigma_1} S_{j \sigma_2} ,

   カノニカル集団の近藤格子模型(Hubbard模型と同様に\ :math:`U`\ と\ :math:`J`\ を入れることも可能)

   .. math::
      :label: fml4_1_kondo

      \mathcal H = - \mu \sum_{i \sigma} c^\dagger_{i \sigma} c_{i \sigma} 
      - t \sum_{\langle i j \rangle \sigma} c^\dagger_{i \sigma} c_{j \sigma} 
      + \frac{J}{2} \sum_{i} \left\{
      S_{i}^{+} c_{i \downarrow}^\dagger c_{i \uparrow}
      + S_{i}^{-} c_{i \uparrow}^\dagger c_{i \downarrow}
      + S_{i z} (n_{i \uparrow} - n_{i \downarrow})
      \right\},

   グランドカノニカル集団のフェルミ粒子Hubbard模型[式 :eq:`fml4_1_hubbard` ]、
   グランドカノニカル集団のスピン模型[式 :eq:`fml4_1_spin` ]、
   グランドカノニカル集団の近藤格子模型[式 :eq:`fml4_1_kondo` ]に対応します。

   ``"SpinGCCMA"``\ では\ ``"SpinGC"``\ と同じ計算を
   より速いアルゴリズム [#]_ を用いて行います。
   ただし、扱うことのできるモデルやMPI並列数に強い制約があります。
   以下の\ ``"Lattice"``\ の項もご参照ください。

-  ``method``

   **形式 :** 文字列(\ ``"Lanczos"``, ``"TPQ"``, ``"Full Diag"``,
   ``"CG"``, ``"Time-Evolution"``\ のいずれか)

   **説明 :** 実行する計算の種類を指定します。
   上記の文字列はそれぞれランチョス法による少数固有状態の計算,
   熱力学的純粋状態を用いた有限温度計算, 直接法による全固有状態計算,
   LOBCG法 [#]_ [#]_ による少数固有状態の計算,
   実時間発展計算 に対応します。

   後述のスペクトル計算において使用される手法もこのパラメーターで指定されます
   ``"CG"``\ とした場合には
   付属している\ :math:`K\omega`\ ライブラリ [#]_ が呼び出され、
   シードスイッチ [#]_ 付きシフト双共役勾配法 [#]_ が適用されます。

*  ``lattice``

   **形式 :** 文字列(\ ``"Chain Lattice"``, ``"Square Lattice"``,
   ``"Triangular Lattice"``, ``"Honeycomb Lattice"``, ``"Ladder"``,
   ``"Kagome"``\ のいずれか)

   **説明 :** 格子の形状を指定します。 上記文字列はそれぞれ1次元鎖(
   :numref:`fig_chap04_1_lattice` (a))、 2次元正方格子(
   :numref:`fig_chap04_1_lattice` (b))、 2次元三角格子(
   :numref:`fig_chap04_1_lattice` (c))、 2次元異方的蜂の巣格子(
   :numref:`fig_chap04_1_honeycomb`)、 梯子格子(:numref:`fig_ladder`)、
   カゴメ格子(:numref:`fig_kagome`)に対応します。

   ``method="SpinGCCMA"``\ では、 このうち\ ``"Chain Lattice"``,
   ``"Honeycomb Lattice"``, ``"Ladder"``,
   ``"Kagome"``\ に対応しています。
   各格子についてのサイズ(\ :math:`L`,\ :math:`W`)とMPI並列数(\ :math:`N_{\rm proc}`)の制限は次のとおりです
   (次節の``L``, ``W``\ もご参照ください)。

   *  ``"Chain Lattice"``

      :math:`L = 8n`\ (ただし:math:`n`\ は\ :math:`n\geq1`\ の整数),
      :math:`N_{\rm proc} \leq 2(L=8)`,
      :math:`N_{\rm proc} \leq 2^{L/2-2}(L>8)`.

   *  ``"Honeycomb Lattice"``

      :math:`W=3, L \geq 2`, :math:`N_{\rm proc} \leq 2(L=2)`,
      :math:`N_{\rm proc} \leq 64(L>2)`.

   *  ``"Ladder"``

      :math:`W=2, L = 2n`\ (ただし:math:`n`\ は\ :math:`n\geq4`\ の整数),
      :math:`N_{\rm proc} \leq 2^{L-4}`.

   *  ``"Kagome"``

      :math:`W=3, L \geq 2`, :math:`N_{\rm proc} \leq 1(L=2)`,
      :math:`N_{\rm proc} \leq 512(L>2)`.

.. [#] \GC=Grand Canonical
.. [#] \Y. Yamaji *et. al.*, manuscript in preparation.
.. [#] A.V.Knyazev, SIAM Journal on Scientific Computing **23**, 517 (2001).
.. [#] S.Yamada, T.Imamura, M.Machida, The Japan Society for Computational Engineering and Science **2006**, 20060027 (2006).
.. [#] https://github.com/issp-center-dev/Komega.
.. [#] S.Yamamoto, T. Sogabe, T. Hoshi, S.-L. Zhang, T. Fujiwara, Journal of the Physical Society of Japan **77**, 114713 (2008).
.. [#] A.Frommer, Computing **70**, 87{109 (2003).


.. raw:: latex

   \newpage
