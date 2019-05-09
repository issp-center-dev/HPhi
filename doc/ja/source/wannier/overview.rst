概要
====

本資料では,
`RESPACK <https://sites.google.com/view/kazuma7k6r>`_ と
mVMC および :math:`{\mathcal H}\Phi` を用いて,
ダウンフォールディングをした格子モデルを計算する機能について説明する.
RESPACKでは, 遮蔽直接積分 :math:`U_{mn}({\bf R},\omega)` および遮蔽交換積分 :math:`J_{mn}({\bf R},\omega)` は以下のような形式で与えられる:

.. math::
   \begin{aligned}
   U_{mn}({\bf R},\omega)&=&\int_V d{\bf r} \int_V  d{\bf r'}
   w_{m{\bf 0}}^*({\bf r}) w_{m{\bf 0}}({\bf r}) 
   W({\bf r,r'},\omega)
   w_{n{\bf R}}^*({\bf r'}) w_{n{\bf R}}({\bf r'}),\nonumber\\
   J_{mn}({\bf R},\omega)&=&\int_V  d{\bf r}\int_V d{\bf r'}
   w_{m{\bf 0}}^*({\bf r}) w_{n{\bf R}}({\bf r}) 
   W({\bf r,r'},\omega) 
   w_{n{\bf R}}^*({\bf r'}) w_{m{\bf 0}}({\bf r'}). 
   \end{aligned}

ここで, :math:`V` は結晶の体積, :math:`w_ {i {\bf R}}({\bf r})` はセル :math:`\bf R` の :math:`i` 番目のワニエ関数, :math:`W({\bf r,r'}, \omega)` は遮蔽クーロン相互作用をそれぞれ表す. 以下, :math:`\omega=0` の成分のみを考慮する. この時、二体相互作用部分のハミルトニアンは以下の形式で与えられる：

.. math::
   \begin{aligned}
   {\cal H}_{\rm int} &= \frac{1}{2}\sum_{\sigma\rho }\sum_{ij}\sum_{nm} \Bigl[ U_{mn}({\bf R}_{ij})c_{im, \sigma}^{\dagger}c_{jn, \rho}^{\dagger}c_{jn, \rho}c_{im, \sigma}\nonumber\\
   &+ J_{mn}({\bf R}_{ij})(c_{im, \sigma}^{\dagger}c_{jn,\rho}^{\dagger}c_{im,\rho}c_{jn,\sigma} + c_{im, \sigma}^{\dagger}c_{im,\rho}^{\dagger}c_{jn,\rho}c_{jn,\sigma}  )\Bigr],
   \end{aligned}

ただし, :math:`{\bf R}_{ij} \equiv {\bf R}_i-{\bf R}_j` とした. ここで, mVMCおよび :math:`{\mathcal H}\Phi` では, :math:`{c_{i, \sigma}^{\dagger}c_{j, \rho}^{\dagger}c_{k, \rho'}c_{l, \sigma'}}` の型の相互作用の入力には対応していないため, 以下のように書き換えたハミルトニアンで定義される：

.. math::
   \begin{aligned}
   {\cal H}_{\rm int} &= \sum_{i,m} U_{mm}({\bf 0})n_{im,\uparrow} n_{im, \downarrow} +\sum_{(i,m)<(j,n)}U_{mn}({\bf R}_{ij})n_{im}n_{jn}\nonumber\\
   & - \sum_{(i,m)<(j,n)}J_{mn}({\bf R}_{ij})(n_{im, \uparrow}n_{jn,\uparrow}+n_{im, \downarrow}n_{jn,\downarrow}) \nonumber\\
   & + \sum_{(i,m)<(j,n)}J_{mn}({\bf R}_{ij})(c_{im, \uparrow}^{\dagger}c_{jn,\downarrow}^{\dagger}c_{im,\downarrow}c_{jn,\uparrow}+{\rm h.c.}) \nonumber\\
   & + \sum_{(i,m)<(j,n)}J_{mn}({\bf R}_{ij}) (c_{im, \uparrow}^{\dagger}c_{im,\downarrow}^{\dagger}c_{jn,\downarrow}c_{jn,\uparrow} + {\rm h.c.} ).
   \end{aligned}


格子モデルは以下のHamiltonianで定義される：

.. math::

   \begin{aligned}
   {\cal H} &=
   \sum_{m,n, i, j,\sigma}
   \left[t_{mn}({\bf R}_{ij}) - t_{mn}^{\rm DC}({\bf R}_{ij})\right] c_{im \sigma}^{\dagger} c_{jn \sigma}
   + {\cal H}_{int},
   \end{aligned}

ただし, :math:`t_{mn}^{\rm DC}({\bf R}_{ij})` は一体項の補正を表し,

.. math::

   \begin{aligned}
   t_{mm}^{\rm DC}({\bf 0}) &\equiv \alpha U_{mm}({\bf 0}) D_{mm}({\bf 0})
   + \sum_{({\bf R}, n) \neq ({\bf 0}, m)} U_{m n} ({\bf R})D_{nn}({\bf 0})\\
   & - (1-\alpha) \sum_{({\bf R}, n) \neq ({\bf 0}, 0)} J_{m n}({\bf R}) D_{nn}({\bf R}),\\
   t_{mn}^{\rm DC}({\bf R}_{ij}) &\equiv \frac{1}{2} J_{mn}({\bf R}_{ij}) \left(D_{nm}({\bf R}_{ji}) + 2 {\rm Re} [D_{nm}({\bf R}_{ji})]\right)\\
   &-\frac{1}{2}  U_{mn}({\bf R}_{ij}) D_{nm}({\bf R}_{ji}),
   \quad ({\bf R}_{ij}, m) \neq ({\bf 0}, n),
   \\
   D_{mn}({\bf R}_{ij}) &\equiv \sum_{\sigma}
   \left\langle c_{im \sigma}^{\dagger} c_{jn \sigma}\right\rangle_{\rm KS},
   \end{aligned}

で与えられる. ここで, :math:`t_{mm}^{\rm DC}({\bf 0})` は化学ポテンシャルの補正項, :math:`t_{mn}^{\rm DC}({\bf R}_{ij})` はFock項に対する補正項を表す. これらは模型を解く際のダブルカウンティングを避けるために導入され, オプションでON/OFFの切り替えが可能になっている. また, 直接積分 :math:`U_{mn}(\bf{R}_{ij})` および交換積分 :math:`J_{mn}({\bf R}_{ij})` をそれぞれ :math:`\lambda_U, \lambda_J` 倍し調節するためのパラメータも用意されている. 
   
要件
----

`QuantumESPRESSO <http://www.quantum-espresso.org/>`_
もしくは
`xTAPP <http://xtapp.cp.is.s.u-tokyo.ac.jp/>`_
を用いてKohn-Sham軌道を用いたのちに,
RESPACKでWannier関数, 誘電関数, 有効相互作用を計算し,
それらを用いて構成した格子モデルを
mVMC もしくは :math:`{\mathcal H}\Phi`
で計算する.
したがってそれらのプログラムが使用可能である必要がある.
