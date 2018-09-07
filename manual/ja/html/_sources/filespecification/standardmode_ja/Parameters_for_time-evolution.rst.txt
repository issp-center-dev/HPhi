実時間発展計算に関するパラメーター
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*  ``dt``

   **形式 :** 正の実数(デフォルトは\ ``0.1``)

   **説明 :** 時間ステップ幅。

*  ``PumpType``

   **形式 :** 文字列(\ ``"Quench"``, ``"Pulse Laser"``, ``"AC Laser"``,
   ``"DC Laser"``\ のいずれか。デフォルトは\ ``"Quench"``)

   **説明 :** 時間依存ハミルトニアンの種類を指定する。
   ``"Quench"``\ では2体演算子
   :math:`U_{\rm quench} \sum_i n_{i \uparrow} n_{i \downarrow}`\ が加えられる。
   ``"Pulse Laser"``, ``"AC Laser"``, ``"DC Laser"``\ では、
   ホッピング項に
   :math:`-\sum_{i j \sigma} t_{i j} \exp[-i{\bf A}(t) \cdot ({\bf R}_i-{\bf R}_j)/(2\pi)] c_{i \sigma} c_{j \sigma}`
   のように位相因子が付く。
   ここで\ :math:`{\bf A}(t)`\ はベクトルポテンシャルであり、
   ``"Pulse Laser"``\ では
   :math:`{\bf A}(t) = {\bf A}_0 \exp[-(t-t_0)^2/(2 t_{\rm dump}^2)] \cos[\omega (t-t_0)]`\ 、
   ``"AC Laser"``\ では
   :math:`{\bf A}(t) = {\bf A}_0 \sin[\omega (t-t_0)]`\ 、
   ``"DC Laser"``\ では :math:`{\bf A}(t) = {\bf A}_0 t`\ となる。

   また、各時刻でのベクトルポテンシャルと電場を図示するためのファイル
   ``potential.dat``\ が出力される。

*  ``Uquench``

   **形式 :** 実数(デフォルトは\ ``0.0``)

   **説明 :** :math:`U_{\rm quench}`

*  ``freq``

   **形式 :** 実数(デフォルトは\ ``0.1``)

   **説明 :** :math:`\omega`

*  ``tshift``

   **形式 :** 実数(デフォルトは\ ``0.0``)

   **説明 :** :math:`t_0`

*  ``tdump``

   **形式 :** 実数(デフォルトは\ ``0.1``)

   **説明 :** :math:`t_{\rm dump}`

*  ``VecPotW``, ``VecPotL``

   **形式 :** 実数(デフォルトはともに\ ``0.0``)

   **説明 :**
   時刻\ :math:`t=t_0`\ でのベクトルポテンシャル\ :math:`{\bf A}_0`\ を
   逆格子のFractional coordinateで指定する。 逆格子ベクトルはFigs.
   :numref:`fig_chap04_1_lattice`, :numref:`fig_chap04_1_honeycomb`,
   :numref:`fig_ladder`, :numref:`fig_kagome`
   に表されている格子ベクトルと対応するものとなる。

.. raw:: latex

   \newpage
