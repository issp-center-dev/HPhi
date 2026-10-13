.. highlight:: none

TMcomponents.dat
~~~~~~~~~~~~~~~~

Lanczos法のリスタートおよび再計算に用いる三重対角行列の要素と初期ノルムを
保存するテキストファイルです。

ファイル名
^^^^^^^^^^

-  ``##_TMComponents.dat``

``##`` はModParaファイルの ``CDataFileHead`` です。

ファイル形式
^^^^^^^^^^^^

-  1行目: ``step``
-  2行目: ``dnorm``
-  続く ``step`` 行: ``alpha beta``
-  基底状態リスタートの場合のみ、最終行:
   ``# nproc P nthreads T``

``alpha`` と ``beta`` はLanczos三重対角行列の要素
:math:`\alpha_i` と :math:`\beta_i` です。末尾のコメントには、対応する
``recalcvec`` の保存に使用したMPIプロセス数とOpenMPスレッド数を記録します。
基底状態リスタートでは、プロセス数またはスレッド数が異なる場合は拒否します。
この行を持たない旧形式ファイルは警告付きで読み込むため、保存時と同じ実行環境を
使用してください。

スペクトル計算の ``CalcSpec`` ファイルは従来のテキスト形式を維持し、
実行環境のコメント行を含みません。

.. raw:: latex

   \newpage
