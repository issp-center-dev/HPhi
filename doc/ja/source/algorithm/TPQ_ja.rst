.. highlight:: none

熱的純粋量子状態による有限温度計算
----------------------------------

杉浦・清水によって、 少数個（サイズが大きい場合はほぼ一つ）の
波動関数から有限温度の物理量を計算する方法が提案されました [1]_ [2]_ 。
その状態は熱的純粋量子状態(TPQ)と呼ばれています。
TPQはハミルトニアンを波動関数に順次作用させて得られるので、
Lanczos法の技術がそのまま使うことができます。
TPQ状態は次のように与えられます。

.. math::

   |\Phi(\beta)\rangle\equiv\exp(-\frac{\beta\hat{\mathcal H}}{2})|\Phi_{\rm rand}\rangle.

ここで :math:`\beta` は逆温度、:math:`|\Phi_{\rm rand}\rangle`  はランダムな初期ベクトルです。
:math:`|\Phi(\beta)\rangle` の期待値として、有限温度の物理量が計算できることが
示されています.

実際の実装について
~~~~~~~~~~~~~~~~~~

mTPQ状態の構成について
^^^^^^^^^^^^^^^^^^^^^^^^^^
ここでは、とくに計算が簡単な, micro canonical TPQ(mTPQ)の
概要を述べます [1]_。

:math:`|\Phi_{\rm rand}\rangle`\ をあるランダムベクトルとします。
これに\ :math:`(l-{\mathcal H}/N_{s})`\ (:math:`l`\ はある定数、\ :math:`N_{s}`\ はサイト数)を\ :math:`k`\ 回作用させた
（規格化された）ベクトルは次のように与えられます。

.. math::

   \begin{aligned}
   |\Phi_{k}\rangle \equiv \frac{(l-{\mathcal H}/N_{s})|\Phi_{k-1}\rangle}{|(l-{\mathcal H}/N_{s})|\Phi_{k-1}\rangle|}.\end{aligned}

この\ :math:`|\Phi_{k}\rangle`\ がmTPQ状態で、このmTPQ状態に対応する逆温度\ :math:`\beta_{k}`\ は
以下のように内部エネルギー\ :math:`u_{k}`\ から求めることができます。

.. math::

   \begin{aligned}
   \beta_{k}\sim \frac{2k/N_{s}}{l-u_{k}},~~
   u_{k} = \langle \Phi_{k}|{\mathcal H}|\Phi_{k}\rangle/N_{s}.\end{aligned}

そして、任意 [3]_ の物理量\ :math:`\hat{A}`\ の\ :math:`\beta_{k}`\ での平均値は

.. math::

   \begin{aligned}
   \langle \hat{A}\rangle_{\beta_{k}} =  \langle \Phi_{k}|\hat{A}|\Phi_{k}\rangle/N_{s}\end{aligned}

となります。 有限系では最初の乱数ベクトルによる誤差がありますので、
いくつか独立な計算を行って、\ :math:`|\psi_{0}\rangle`
に関する平均値および標準偏差を見積もっています。

cTPQ状態の構成について
^^^^^^^^^^^^^^^^^^^^^^^^^^
カノニカルTPQ(cTPQ)状態の構成方法について述べます [2]_,
cTPQ法では :math:`\exp[-\beta\hat{\mathcal H}/2]` は
次にように近似されます.

.. math::

  &\exp(-\frac{\beta\hat{\mathcal H}}{2})|\Phi_{\rm rand}\rangle\sim|\Phi_{k}\rangle = U(\Delta\tau)^{k}|\Phi_{\rm rand}\rangle\\
  &U(\Delta\tau)=\sum_{n=0}^{n_{\rm max}}\frac{1}{n!}(-\frac{\Delta\tau}{2}\hat{\mathcal H})^{n}\\
  &\beta_{k}=k\Delta \tau
   
mTPQと同じように物理量はcTPQ状態の期待値として計算できます。

.. math::

   \langle \hat{A}\rangle_{\beta_{k}} =  \frac{\langle\Phi_{k}|\hat{A}|\Phi_{k}\rangle}{\langle\Phi_{k}|\Phi_{k}\rangle}.



初期ベクトルの設定について
^^^^^^^^^^^^^^^^^^^^^^^^^^

熱的純粋量子状態による有限温度計算では、初期ベクトルは全ての成分に対してランダムな係数を与えます。
初期ベクトルの係数の型はModParaで指定される入力ファイルの\ ``InitialVecType``\ を用い、
実数もしくは複素数の指定をすることができます。mTPQとcTPQはともに
SIMD-oriented Fast Mersenne Twister (dSFMT) [4]_ を使います。
``initial_iv`` が基準seed、``NumAve`` がsample数、``InitialVecType`` が成分の分布を
指定します。Standard modeの入力、またはExpert modeのModPara/CalcModで指定します。

現在の乱数stream versionは2です。sample/rank/threadごとに、次の7個の32-bit語を
配列としてdSFMTへ渡して初期化します::

   [0x48505451, 2, low32(abs(initial_iv)), high32(abs(initial_iv)),
    sample_index, MPI_rank, OpenMP_thread]

番号は0始まりです。先頭のtagとversionで方式を識別し、sample・rank・threadは単一seedへの
加算値ではなく独立した語として渡します。``initial_iv`` の絶対値は符号なし整数として
全bitを保持し、最小の負整数も扱います。絶対値が等しい正負seedは同じ意味です。
0を指定した場合も全sampleで同じstreamを再利用しません。OpenMPはstatic scheduleで
連続した成分区間を各threadへ割り当てます。

``NumAve`` はstream keyに含めないため、入力と並列構成を固定してsample数を増やしても
既存sampleは変わりません。MPI/OpenMP数の変更は成分の所有・割当を変えるので、初期vectorが
変わる場合があります。浮動小数点による正規化もplatform間で差があり得ます。
新規初期化は実行logに ``TPQ initial RNG: stream_version=2 sample=...`` と記録します。

旧方式は ``123432 + (sample+1)*abs(initial_iv) + thread + threads*rank`` という
加算seedでした。既定の ``abs(initial_iv)=1`` では隣接sampleが別workerの乱数列を
再利用し、0ではsample間で同じ列を再利用していました。version 2はこのseed重複を解消するため、
同じ入力seedでも新規計算の結果は旧版と変わります。保存stateからの継続は乱数を生成し直さず、
引き続き可能です。ただし旧checkpointの継続によって既存sampleの相関は解消されません。
追加する新sampleはversion 2を使います。旧計算の熱平均・誤差棒への影響量は計算条件に依存し、
乱数列重複の回帰試験だけからは確定しません。

.. [1] \S. Sugiura, A. Shimizu, Phys. Rev. Lett. **108**, 240401 (2012).
.. [2] \S. Sugiura and A. Shimizu, Phys. Rev. Lett. **111**, 010401 (2013).
.. [3] 局所的にという条件がつきます。
.. [4] \http://www.math.sci.hiroshima-u.ac.jp/m-mat/MT/SFMT.
