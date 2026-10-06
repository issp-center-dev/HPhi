.. highlight:: none

Time\_TPQ\_Step.dat
~~~~~~~~~~~~~~~~~~~

| (TPQ法でのみ出力) TPQ法でのステップ毎の開始時間を出力します。
  再計算の場合は値が追記されます。 以下にファイル例を示します。

::

    set 0 step 1:TPQ begins: Wed Jul 13 07:59:20 2016
    set 0 step 2:TPQ begins: Wed Jul 13 07:59:20 2016
    set 0 step 3:TPQ begins: Wed Jul 13 07:59:20 2016
    …
    set 4 step 1997:TPQ begins: Wed Jul 13 07:59:32 2016
    set 4 step 1998:TPQ begins: Wed Jul 13 07:59:32 2016
    set 4 step 1999:TPQ begins: Wed Jul 13 07:59:32 2016

ファイル名
^^^^^^^^^^

-  ##\_TPQ\_Step.dat

##はModParaファイル内の[string02]で指定されるヘッダを表します。

ファイル形式
^^^^^^^^^^^^

-  set :math:`[`\ int01\ :math:`]` stp :math:`[`\ int02\ :math:`]`: TPQ
   begins: :math:`[`\ string01\ :math:`]`

パラメータ
^^^^^^^^^^

-  :math:`[`\ int01\ :math:`]`

   **形式 :** int型

   **説明 :** TPQ法でのシード数。

-  :math:`[`\ int02\ :math:`]`

   **形式 :** int型

   **説明 :** TPQ法でのステップ数。

-  :math:`[`\ string01\ :math:`]`

   **形式 :** string型

   **説明 :** 計算開始時間。

.. raw:: latex

   \newpage