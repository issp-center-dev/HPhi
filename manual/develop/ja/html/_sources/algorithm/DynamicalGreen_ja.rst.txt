.. highlight:: none

動的グリーン関数
----------------

:math:`{\cal H}\Phi`\ では励起状態\ :math:`|\Phi ' \rangle  = \hat{O} | \Phi _0 \rangle`\ に対する動的関数

.. math:: I(z) = \langle \Phi ' | \frac{1}{ {\cal H}- z\hat{I} } | \Phi '\rangle

を計算することができます。 演算子\ :math:`\hat{O}`\ はシングル励起状態

.. math:: \sum_{i, \sigma_1} A_{i \sigma_1} c_{i \sigma_1} (c_{i\sigma_1}^{\dagger})

およびペア励起状態

.. math:: \sum_{i, j, \sigma_1, \sigma_2} A_{i \sigma_1 j \sigma_2} c_{i \sigma_1}c_{j \sigma_2}^{\dagger} (c_{i\sigma_1}^{\dagger}c_{j\sigma_2})


として、それぞれ定義することが出来ます。例えば、動的スピン感受率を計算する場合はペア励起演算子を用い

.. math:: \hat{O} = \hat{S}({\bf k}) = \sum_{j}\hat{S}_j^z e^{i  {\bf k} \cdot \bf {r}_j} = \sum_{j}\frac{1}{2} (c_{j\uparrow}^{\dagger}c_{j\uparrow}-c_{j\downarrow}^{\dagger}c_{j\downarrow})e^{i  {\bf k} \cdot \bf {r}_j}

のように定義することで計算することができます。
なお、動的グリーン関数の計算には、Lanczos法を用いた連分数展開による解法 [1]_ とシフト型クリロフ理論による解法 [2]_ の2つが実装されています。
詳細については各文献を参照してください。

.. [1] \E. daggerotto, Rev. Mod. Phys. **66**, 763-840 (1994).
.. [2] \S.Yamamoto, T. Sogabe, T. Hoshi, S.-L. Zhang, T. Fujiwara, Journal of the Physical Society of Japan **77**, 114713 (2008).