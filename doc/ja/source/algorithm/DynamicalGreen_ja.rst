.. highlight:: none

動的Green関数
-------------

:math:`{\cal H}\Phi`\ では動的関数

.. math:: G_n^{O_l,O_r}(z) = \langle \Phi_n | \hat{O}_l (z + E_n - \hat{\cal H})^{-1} \hat{O}_r| \Phi_n \rangle

を計算することができます。 演算子\ :math:`\hat{O}_{l,r}`\ はシングル励起状態

.. math:: \sum_{i, \sigma_1} A_{i \sigma_1} c_{i \sigma_1} \quad \textrm{or} \quad \sum_{i, \sigma_1} A_{i \sigma_1} c_{i\sigma_1}^{\dagger}

およびペア励起状態

.. math:: \sum_{i, j, \sigma_1, \sigma_2} A_{i \sigma_1 j \sigma_2} c_{i \sigma_1}c_{j \sigma_2}^{\dagger} \quad \textrm{or} \quad
          \sum_{i, j, \sigma_1, \sigma_2} A_{i \sigma_1 j \sigma_2} c_{i\sigma_1}^{\dagger}c_{j\sigma_2}


として、それぞれ定義することが出来ます。例えば、動的スピン感受率を計算する場合はペア励起演算子を用い

.. math:: \hat{O}_r = \hat{S}_{\textbf{R}=\textbf{0}}^z = \frac{1}{2} (c_{\textbf{0}\uparrow}^{\dagger}c_{\textbf{0}\uparrow}-c_{\textbf{0}\downarrow}^{\dagger}c_{\textbf{0}\downarrow})
    \\
    \hat{O}_l = \hat{S}_{\textbf{R}}^z = \frac{1}{2} (c_{\textbf{R}\uparrow}^{\dagger}c_{\textbf{R}\uparrow}-c_{\textbf{R}\downarrow}^{\dagger}c_{\textbf{R}\downarrow})

として、:math:`G_n^{O_l,O_r}(z)\equiv G_n^{\textbf{R}}(z)` を計算し、ポストプロセスで

.. math:: G_n^{\textbf{k}}(z) \equiv \sum_{\textbf{R}} \exp(i\textbf{k}\cdot\textbf{R}) G_n^{\textbf{R}}(z)

のようにFourier変換を行い計算することができます。
なお、動的関数の計算には、Lanczos法を用いた連分数展開による解法 [1]_ 、シフト型クリロフ理論による解法 [2]_ 、およびLehmann表示による動的関数

.. math:: G_n^{O_l,O_r}(z) = \sum_{m} \frac{\langle \Phi_n | \hat{O}_l | \Phi_m \rangle \langle \Phi_m |\hat{O}_r| \Phi_n \rangle}{z + E_n - E_m}

を全対角化により直接計算する手法の3つが実装されています。
詳細については各文献を参照してください。

.. [1] \E. Dagotto, Rev. Mod. Phys. **66**, 763-840 (1994).
.. [2] \S.Yamamoto, T. Sogabe, T. Hoshi, S.-L. Zhang, T. Fujiwara, Journal of the Physical Society of Japan **77**, 114713 (2008).
