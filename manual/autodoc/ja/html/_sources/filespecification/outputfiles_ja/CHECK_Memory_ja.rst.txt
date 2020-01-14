.. highlight:: none

CHECK\_Memory.dat
~~~~~~~~~~~~~~~~~

使用するメモリサイズの出力を行います。配列サイズおよび必要なメモリを出力します。
以下にファイル例を記載します。

::

    MAX DIMENSION idim_max=400 
    REQUIRED MEMORY  max_mem=0.000019 GB 

ファイル形式
^^^^^^^^^^^^

以下のようなファイル形式をとります。

-  MAX DIMENSION idim\_max=\ :math:`[`\ int01\ :math:`]`

-  REQUIRED MEMORY max\_mem =\ :math:`[`\ double01\ :math:`]` GB

パラメータ
^^^^^^^^^^

-  :math:`[`\ int01\ :math:`]`

   **形式 :** int型

   **説明 :** ヒルベルトスペースの数を表します。

-  :math:`[`\ double01\ :math:`]`

   **形式 :** double型

   **説明 :**
   ヒルベルトスペースの確保に必要とするメモリサイズを表します(単位はGB)。

.. raw:: latex

   \newpage