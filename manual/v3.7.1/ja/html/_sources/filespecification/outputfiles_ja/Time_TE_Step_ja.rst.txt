.. highlight:: none

Time\_TE\_Step.dat
~~~~~~~~~~~~~~~~~~

| (実時間発展法でのみ出力)
  実時間発展法でのステップ毎の開始時間を出力します。以下にファイル例を示します。

::

    step 1:TE begins: Wed Jul 13 07:59:20 2017
    step 2:TE begins: Wed Jul 13 07:59:20 2017
    step 3:TE begins: Wed Jul 13 07:59:20 2017
    …
    step 1997:TE begins: Wed Jul 13 07:59:32 2017
    step 1998:TE begins: Wed Jul 13 07:59:32 2017
    step 1999:TE begins: Wed Jul 13 07:59:32 2017

ファイル名
^^^^^^^^^^

-  ##\_TE\_Step.dat

##はModParaファイル内の[string02]で指定されるヘッダを表します。

ファイル形式
^^^^^^^^^^^^

-  stp :math:`[`\ int01\ :math:`]`: TE begins:
   :math:`[`\ string01\ :math:`]`

パラメータ
^^^^^^^^^^

-  :math:`[`\ int01\ :math:`]`

   **形式 :** int型

   **説明 :** TE法でのステップ数。

-  :math:`[`\ string01\ :math:`]`

   **形式 :** string型

   **説明 :** 計算開始時間。


.. raw:: latex

   \newpage