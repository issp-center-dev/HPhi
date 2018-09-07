.. highlight:: none

.. _Subsec:InputFileList:


入力ファイル指定用ファイル
~~~~~~~~~~~~~~~~~~~~~~~~~~

| 計算で使用する入力ファイル一式を指定します。ファイル形式に関しては、以下のようなフォーマットをしています。

::

    CalcMod  calcmod.def
    ModPara  modpara.def
    LocSpin  zlocspn.def
    Trans    ztransfer.def
    InterAll zinterall.def
    OneBodyG zcisajs.def
    TwoBodyG    zcisajscktaltdc.def

| 

ファイル形式
^^^^^^^^^^^^

[string01] [string02]

パラメータ
^^^^^^^^^^

-  :math:`[`\ string01\ :math:`]`

   **形式 :** string型 (固定)

   **説明 :** キーワードを指定します。

-  :math:`[`\ string02\ :math:`]`

   **形式 :** string型

   **説明 :** キーワードにひも付けられるファイル名を指定します(任意)。

使用ルール
^^^^^^^^^^

本ファイルを使用するにあたってのルールは以下の通りです。

-  キーワードを記載後、半角空白を開けた後にファイル名を書きます。ファイル名は自由に設定できます。

-  必ず指定しなければいけないパラメーターはCalcMod, ModPara,
   LocSpinです。

-  各キーワードは順不同に記述できます。

-  指定したキーワード、ファイルが存在しない場合はエラー終了します。

-  :math:`\#`\ で始まる行は読み飛ばされます。

-  ファイル読込用キーワードは :numref:`Table 4.2` により指定します。

.. _Table 4.2:
.. csv-table:: 定義ファイル一覧
    :header: "Keywords", "指定ファイルの概要"
    :widths: 4, 20

    "CalcMod","計算モードに関する指定をします。"
    "ModPara","計算で用いるパラメータの指定をします。"
    "LocSpin","各サイトに対して遍歴電子もしくは局在スピンの指定をします。"
    "Trans","一般的一体相互作用に関する設定をします。"
    "InterAll", "一般的二体相互作用に関する設定をします。"
    "CoulombIntra", "内部クーロン相互作用に関する設定をします。"
    "CoulombInter", "サイト間クーロン相互作用に関する設定をします。"
    "Hund", "フント結合に関する設定をします。"
    "PairHop", "ペアホッピングに関する設定をします。"
    "Exchange", "交換相互作用に関する設定をします。"
    "Ising", "イジング相互作用に関する設定をします。"
    "PairLift", "ペアリフト相互作用に関する設定をします。"
    "OneBodyG", "出力する一体グリーン関数 \ :math:`\langle c_{i\sigma}^{\dagger}c_{j\sigma}\rangle` に関する設定をします。"
    "TwoBodyG", "出力するニ体グリーン関数 :math:`\langle c_{i\sigma}^{\dagger}c_{j\sigma}c_{k\tau}^{\dagger}c_{l\tau}\rangle`\ に関する設定をします。"
    "SingleExcitation", "一体励起状態の生成演算子に関する指定をします。"
    "PairExcitation", "ニ体励起状態の生成演算子に関する指定をします。"
    "SpectrumVec", "スペクトル関数を計算するためのリスタート用の入力ベクトルを指定します。"                                                               

.. raw:: latex

   \newpage
