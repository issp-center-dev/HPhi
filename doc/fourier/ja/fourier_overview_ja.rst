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
   \langle {\hat c}_{{\bf k} \uparrow}^{\dagger} {\hat c}_{{\bf k} \uparrow}\rangle
   &\equiv \frac{1}{N_{\rm cell}} \sum_{i j}^{N_{\rm site}} e^{-i {\bf k}\cdot({\bf R}_i - {\bf R}_j)}
   \langle {\hat c}_{i \uparrow}^{\dagger} {\hat c}_{j \uparrow}\rangle
   \\
   \langle {\hat c}_{{\bf k} \downarrow}^{\dagger} {\hat c}_{{\bf k} \downarrow}\rangle
   &\equiv \frac{1}{N_{\rm cell}} \sum_{i j}^{N_{\rm site}} e^{-i {\bf k}\cdot({\bf R}_i - {\bf R}_j)}
   \langle {\hat c}_{i \downarrow}^{\dagger} {\hat c}_{j \downarrow}\rangle
   \end{align}

密度-密度相関

.. math::

   \begin{align}
   \langle {\hat \rho}_{\bf k} {\hat \rho}_{\bf k}\rangle
   \equiv \frac{1}{N_{\rm cell}^2} \sum_{i j}^{N_{\rm site}} e^{-i {\bf k}\cdot({\bf R}_i - {\bf R}_j)}
   \langle ({\hat \rho}_{i} - \langle {\hat \rho}_{i} \rangle)
           ({\hat \rho}_{j} - \langle {\hat \rho}_{j} \rangle) \rangle
   \end{align}

スピン-スピン相関

.. math::
   :nowrap:

   \begin{align}
   \langle {\hat S}_{\bf k}^{z} {\hat S}_{\bf k}^{z} \rangle
   &\equiv \frac{1}{N_{\rm cell}^2} \sum_{i j}^{N_{\rm site}} e^{-i {\bf k}\cdot({\bf R}_i - {\bf R}_j)}
   \langle {\hat S}_{i}^{z} {\hat S}_{j}^{z} \rangle
   \\
   \langle {\hat S}_{\bf k}^{+} {\hat S}_{\bf k}^{-} \rangle
   &\equiv \frac{1}{N_{\rm cell}^2} \sum_{i j}^{N_{\rm site}} e^{-i {\bf k}\cdot({\bf R}_i - {\bf R}_j)}
   \langle {\hat S}_{i}^{+} {\hat S}_{j}^{-} \rangle
   \\
   \langle {\hat {\bf S}}_{\bf k} \cdot {\hat {\bf S}}_{\bf k} \rangle
   &\equiv \frac{1}{N_{\rm cell}^2} \sum_{i j}^{N_{\rm site}} e^{-i {\bf k}\cdot({\bf R}_i - {\bf R}_j)}
   \langle {\hat {\bf S}}_{i} \cdot {\hat {\bf S}}_{j} \rangle
   \end{align}

