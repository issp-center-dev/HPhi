.. highlight:: none

WarningOnTransfer.dat
~~~~~~~~~~~~~~~~~~~~~

Transferの成分が重複している場合に出力されます。
以下にファイル例を記載します。

::

    double conuntings in transfers: i=0 j=2 spni 0 spnj 0  
    double conuntings in transfers: i=2 j=0 spni 0 spnj 0  
    double conuntings in transfers: i=0 j=2 spni 1 spnj 1  
    double conuntings in transfers: i=2 j=0 spni 1 spnj 1  

ファイル形式
^^^^^^^^^^^^

以下のようなファイル形式をとります。

-  double conuntings in transfers: i=\ :math:`[`\ int01\ :math:`]`
   j=\ :math:`[`\ int02\ :math:`]` spni :math:`[`\ int03\ :math:`]` spnj
   :math:`[`\ int04\ :math:`]`

パラメータ
^^^^^^^^^^

-  :math:`[`\ int01\ :math:`]`, :math:`[`\ int02\ :math:`]`

   **形式 :** int型

   **説明 :** 重複しているTransferのサイト番号を表します。

-  :math:`[`\ int03\ :math:`]`, :math:`[`\ int04\ :math:`]`

   **形式 :** int型

   | **説明 :** 重複しているTransferのスピン番号を表します。
   | 0: アップスピン
   | 1: ダウンスピン
   | を表します。

.. raw:: latex

   \newpage
