.. highlight:: none

実時間発展
----------

:math:`{\cal H}\Phi`\ では初期波動関数を\ :math:`|\Phi(t_0)\rangle`\ として、

.. math:: |\Phi (t_n)\rangle = e^{-i {\cal H}  \Delta t_n}|\Phi (t_{n-1})\rangle


の関係を用いて逐次時間発展の計算をしています。
ここで時刻 :math:`t_n = \sum_{j=1}^n  \Delta t_j` です。
実際の計算では\ :math:`e^{-i {\cal H}  \Delta t_n}`\ を

.. math:: e^{-i {\cal H}  \Delta t_n} =\sum_{l=0}^m \frac{1}{l!}(-i {\cal H}  \Delta t_n)^l


と近似しています。展開次数の上限\ :math:`m`\ は\ ``ModPara``\ ファイルの\ ``ExpandCoef``\ で指定することができます。
展開次数が十分かどうかは、
ノルム\ :math:`\langle \Phi (t_n)|\Phi (t_n)\rangle=1`\ が成立しているかで検証することができます。