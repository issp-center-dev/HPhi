概要
====

本資料は, mVMC および :math:`{\mathcal H}\Phi` で計算された
サイト表示の相関関数をFourier変換し, 出力するユーティリティに関するマニュアルである.

要件
----

本ユーティリティの使用要件はmVMC および :math:`{\mathcal H}\Phi` と同じである.

.. _supported:

対応する量
----------

本ユーティリティは以下の相関関数のFourier変換に対応している.

1体相関

.. math::
   :nowrap:

   \begin{align}
   \langle {\hat c}_{{\bf k} \alpha \uparrow}^{\dagger} {\hat c}_{{\bf k} \beta \uparrow}\rangle
   &\equiv \sum_{\bf R}^{N_{\bf R}} e^{-i {\bf k}\cdot{\bf R}}
   \langle {\hat c}_{{\bf 0} \alpha \uparrow}^{\dagger} {\hat c}_{{\bf R} \beta \uparrow}\rangle
   \\
   \langle {\hat c}_{{\bf k} \alpha \downarrow}^{\dagger} {\hat c}_{{\bf k} \beta \downarrow}\rangle
   &\equiv \sum_{\bf R}^{N_{\bf R}} e^{-i {\bf k}\cdot {\bf R}}
   \langle {\hat c}_{{\bf 0} \alpha \downarrow}^{\dagger} {\hat c}_{{\bf R} \beta \downarrow}\rangle
   \end{align}

密度-密度相関

.. math::

   \begin{align}
   \langle {\hat \rho}_{{\bf k}\alpha} {\hat \rho}_{{\bf k}\beta}\rangle
   \equiv \frac{1}{N_{\bf R}} \sum_{\bf R}^{N_{\bf R}} e^{-i {\bf k}\cdot{\bf R}}
   \langle ({\hat \rho}_{{\bf 0}\alpha} - \langle {\hat \rho}_{{\bf 0}\alpha} \rangle)
           ({\hat \rho}_{{\bf R}\beta} - \langle {\hat \rho}_{{\bf R}\beta} \rangle) \rangle
   \end{align}

スピン-スピン相関

.. math::
   :nowrap:

   \begin{align}
   \langle {\hat S}_{{\bf k}\alpha}^{z} {\hat S}_{{\bf k}\beta}^{z} \rangle
   &\equiv \frac{1}{N_{\bf R}} \sum_{\bf R}^{N_{\bf R}} e^{-i {\bf k}\cdot{\bf R}}
   \langle {\hat S}_{{\bf 0}\alpha}^{z} {\hat S}_{{\bf R}\beta}^{z} \rangle
   \\
   \langle {\hat S}_{{\bf k}\alpha}^{+} {\hat S}_{{\bf k}\beta}^{-} \rangle
   &\equiv \frac{1}{N_{\bf R}} \sum_{\bf R}^{N_{\bf R}} e^{-i {\bf k}\cdot{\bf R}}
   \langle {\hat S}_{{\bf 0}\alpha}^{+} {\hat S}_{{\bf R}\beta}^{-} \rangle
   \\
   \langle {\hat {\bf S}}_{{\bf k}\alpha} \cdot {\hat {\bf S}}_{{\bf k}\beta} \rangle
   &\equiv \frac{1}{N_{\bf R}} \sum_{\bf R}^{N_{\bf R}} e^{-i {\bf k}\cdot{\bf R}}
   \langle {\hat {\bf S}}_{{\bf 0}\alpha} \cdot {\hat {\bf S}}_{{\bf R}\beta} \rangle
   \end{align}
