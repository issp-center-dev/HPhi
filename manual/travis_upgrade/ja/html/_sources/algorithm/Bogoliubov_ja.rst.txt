.. highlight:: none

.. _Sec:sec_bogoliubov_rep:

Bogoliubov表現
--------------

スピン系の計算において一体項(\ ``transfer``)、\ ``InterAll``\ 形式での相互作用、
相関関数のインデックスの指定にはBogoliubov表現が使われています。
スピンの演算子は次のように生成\ :math:`\cdot`\ 消滅演算子で書き換えられます。

.. math::

  S_{i z} &= \sum_{\sigma = -S}^{S} \sigma c_{i \sigma}^\dagger c_{i \sigma}
  \\
  S_{i}^+ &= \sum_{\sigma = -S}^{S-1} 
  \sqrt{S(S+1) - \sigma(\sigma+1)} 
  c_{i \sigma+1}^\dagger c_{i \sigma}
  \\
  S_{i}^- &= \sum_{\sigma = -S}^{S-1} 
  \sqrt{S(S+1) - \sigma(\sigma+1)} 
  c_{i \sigma}^\dagger c_{i \sigma+1}


