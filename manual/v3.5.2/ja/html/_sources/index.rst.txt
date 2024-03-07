#####################################################
:math:`{\mathcal H}\Phi` のドキュメントへようこそ！
#####################################################

:math:`{\mathcal H}\Phi` とは？
------------------------------------------
並列計算機に対応した数値厳密対角化法による有効模型ソルバーパッケージ。広汎な多体量子系の有効模型(多軌道ハバード模型、ハイゼンベルグ模型、近藤格子模型など)の基底状態及び低励起状態の波動関数を並列計算によって求められます。ランチョス法による基底状態計算、熱的純粋量子状態を利用した比熱・帯磁率の温度依存性計算が可能です。さらに、ver.2.0より数値ライブラリKωが接続され、先端的数理アルゴリズム(シフト型クリロフ部分空間理論)による動的グリーン関数の計算が可能となっています。

ライセンス
--------------
本ソフトウェアのプログラムパッケージおよびソースコード一式はGNU General Public License version 3（GPL v3）に準じて配布されています。

:math:`{\mathcal H}\Phi` (hphi)を引用する際には、以下の文献を引用してください。

Mitsuaki Kawamura, Kazuyoshi Yoshimi, Takahiro Misawa, Youhei Yamaji, Synge Todo, and Naoki Kawashima

    `Comp. Phys. Commun. 217 (2017) 180-192 <http://www.sciencedirect.com/science/article/pii/S0010465517301200>`_.


コピーライト
------------------

© *2015- The University of Tokyo. All rights reserved.*

本ソフトウェアは2015, 2016, 2017年度 東京大学物性研究所 ソフトウェア高度化プロジェクトの支援を受け開発されており、その著作権は東京大学が所持しています。

ダウンロード
------------------
:math:`{\mathcal H}\Phi` のソースコードは `GitHub page <https://github.com/issp-center-dev/HPhi>`_ or `release page <https://github.com/issp-center-dev/HPhi/releases>`_ からダウンロードできます。


Contents
--------
.. toctree::
   :maxdepth: 3
   :numbered: 3

   introduction_ja
   howtouse/ho-index
   tutorial/tu-index
   filespecification/fi-index
   algorithm/al-index
   acknowledgement_ja
   fourier/index
   wannier/index

.. Indices and tables
.. ==================

.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`
