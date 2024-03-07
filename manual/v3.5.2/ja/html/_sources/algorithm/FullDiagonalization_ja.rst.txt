.. highlight:: none

完全対角化
----------

手法概要
~~~~~~~~

ハミルトニアン\ :math:`{{\mathcal H}}`\ を実空間配置\ :math:`| \psi_j \rangle`\ (:math:`j=1\cdots N`)を用いて
作成します:
:math:`H_{ij}= \langle \psi_i | {\hat H} | \psi_j \rangle`\ 。
この行列を対角化することで、固有値\ :math:`E_i`\ 、
固有ベクトル\ :math:`|\Phi_i\rangle`\ を求める
ことができます(\ :math:`i=1 \cdots N`)。なお、対角化ではlapackの\ **dsyev**\ また
**zheev**\ を用いています。 また、有限温度計算用に各固有エネルギー状態の
期待値\ :math:`\langle A_i\rangle \equiv \langle \Phi_i | {\hat A} | \Phi_i\rangle`\ を計算・出力するようにしています。

有限温度物理量の計算
~~~~~~~~~~~~~~~~~~~~

完全対角化で求めた\ :math:`\langle A_i\rangle \equiv \langle \Phi_i | {\hat A} | \Phi_i\rangle`\ を用い、

.. math:: \langle {\hat A}\rangle=\frac{\sum_{i=1}^N \langle A_i\rangle {\rm  e}^{-\beta E_i}}{\sum_{i=1}^N{\rm  e}^{-\beta E_i}}

の関係から有限温度の物理量を計算します。
実際の計算処理としてはポスト処理により計算を行います。