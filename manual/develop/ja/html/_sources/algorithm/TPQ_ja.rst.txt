.. highlight:: none

熱的純粋量子状態による有限温度計算
----------------------------------

杉浦・清水によって、 少数個（サイズが大きい場合はほぼ一つ）の
波動関数から有限温度の物理量を計算する方法が提案されました [1]_ 。
その状態は熱的純粋量子状態(TPQ)と呼ばれています。
TPQはハミルトニアンを波動関数に順次作用させて得られるので、
Lanczos法の技術がそのまま使うことができます。
ここでは、とくに計算が簡単な, micro canonical TPQ(mTPQ)の
概要を述べます。

:math:`|\psi_{0}\rangle`\ をあるランダムベクトルとします。
これに\ :math:`(l-{\mathcal H}/N_{s})`\ (:math:`l`\ はある定数、\ :math:`N_{s}`\ はサイト数)を\ :math:`k`\ 回作用させた
（規格化された）ベクトルは次のように与えられます。

.. math::

   \begin{aligned}
   |\psi_{k}\rangle \equiv \frac{(l-{\mathcal H}/N_{s})|\psi_{k-1}\rangle}{|(l-{\mathcal H}/N_{s})|\psi_{k-1}\rangle|}.\end{aligned}

この\ :math:`|\psi_{k}\rangle`\ がmTPQ状態で、このmTPQ状態に対応する逆温度\ :math:`\beta_{k}`\ は
以下のように内部エネルギー\ :math:`u_{k}`\ から求めることができます。

.. math::

   \begin{aligned}
   \beta_{k}\sim \frac{2k/N_{s}}{l-u_{k}},~~
   u_{k} = \langle \psi_{k}|{\mathcal H}|\psi_{k}\rangle/N_{s}.\end{aligned}

そして、任意 [2]_ の物理量\ :math:`\hat{A}`\ の\ :math:`\beta_{k}`\ での平均値は

.. math::

   \begin{aligned}
   \langle \hat{A}\rangle_{\beta_{k}} =  \langle \psi_{k}|\hat{A}|\psi_{k}\rangle/N_{s}\end{aligned}

となります。 有限系では最初の乱数ベクトルによる誤差がありますので、
いくつか独立な計算を行って、\ :math:`|\psi_{0}\rangle`
に関する平均値および標準偏差を見積もっています。

実際の実装について
~~~~~~~~~~~~~~~~~~

初期ベクトルの設定について
^^^^^^^^^^^^^^^^^^^^^^^^^^

熱的純粋量子状態による有限温度計算では、初期ベクトルは全ての成分に対してランダムな係数を与えます。
初期ベクトルの係数の型はModParaで指定される入力ファイルの\ ``InitialVecType``\ を用い、
実数もしくは複素数の指定をすることができます。乱数のシードは\ ``initial_iv``
(:math:`\equiv r_s`)により

.. math::

   \begin{aligned}
   123432+(n_{\rm run}+1)\times  |r_s|+k_{\rm Thread}+N_{\rm Thread} \times k_{\rm Process}\end{aligned}


で与えられます。ここで、\ :math:`n_{\rm run}`\ はrunの回数であり、runの最大回数はスタンダードモード用入力ファイル、
もしくはModParaで指定される入力ファイルの\ ``NumAve``\ で指定します。
``initial_iv``\ はスタンダードモード用の入力ファイル、もしくはエキスパートモードではModParaで指定される入力ファイルで指定します。乱数はSIMD-oriented
Fast Mersenne Twister(dSFMT)を用い発生させています [3]_ 。
また、\ :math:`k_{\rm Thread}, N_{\rm Thread}, k_{\rm Process}`\ 
はそれぞれスレッド番号、スレッド数、プロセス番号を表します。
したがって同じ\ ``initial_iv``\ を用いても、並列数が異なる場合には別の初期波動関数が生成されます。

.. [1] \S. Sugiura, A. Shimizu, Phys. Rev. Lett. **108**, 240401 (2012).
.. [2] 局所的にという条件がつきます。
.. [3] \http://www.math.sci.hiroshima-u.ac.jp/m-mat/MT/SFMT.