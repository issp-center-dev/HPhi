概要
====

本資料では,
`RESPACK <https://sites.google.com/view/kazuma7k6r>`_ と
mVMC および :math:`{\mathcal H}\Phi` を用いて,
ダウンフォールディングをした格子モデルを計算する機能について説明する.

.. math::

   \begin{aligned}
   {\cal H} &=
   \sum_{i, j, \alpha, \beta, \sigma}
   t_{i \alpha j \beta} c_{i \alpha \sigma}^{\dagger} c_{j \beta \sigma}
   \nonumber \\
   &+ \sum_{i, \alpha}
   U_{i \alpha i \alpha} n_{i \alpha \uparrow} n_{j \alpha \downarrow}
   + \sum_{(i, \alpha) < (j, \beta)}
   U_{i \alpha j \beta} n_{i \alpha} n_{j \beta}
   - \sum_{(i, \alpha) < (j, \beta)}
   J_{i \alpha j \beta} (n_{i \alpha \uparrow} n_{j \beta \uparrow}
   + n_{i \alpha \downarrow} n_{j \beta \downarrow})
   \nonumber \\
   &+ \sum_{(i, \alpha) < (j, \beta)}
   J_{i \alpha j \beta} (
   c_{i \alpha \uparrow}^{\dagger} c_{j \beta \downarrow}^{\dagger}
   c_{i \alpha \downarrow} c_{j \beta \uparrow} +
   c_{j \beta \uparrow}^{\dagger} c_{i \alpha \downarrow}^{\dagger}
   c_{j \beta \downarrow} c_{j \alpha \uparrow} )
   \nonumber \\
   &+ \sum_{(i, \alpha) < (j, \beta)}
   J_{i \alpha j \beta} (
   c_{i \alpha \uparrow}^{\dagger} c_{i \alpha \downarrow}^{\dagger}
   c_{j \beta \downarrow} c_{j \beta \uparrow} +
   c_{j \beta \uparrow}^{\dagger} c_{j \beta \downarrow}^{\dagger}
   c_{i \alpha \downarrow} c_{i \alpha \uparrow} )
   \end{aligned}

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
