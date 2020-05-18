.. _Ch:HowToStandard:

.. highlight:: none

スタンダードモード用入力ファイル
================================

スタンダードモード用入力ファイルは次のような格好をしています。
::

    W = 2
    L = 4
    model = "spin"
    method = "Lanczos"

    lattice = "triangular lattice"
    //mu = 1.0
    // t = -1.0
    // t' = -0.5
    // U = 8.0
    //V = 4.0
    //V'=2.0
    J = -1.0
    J'=-0.5
    // nelec = 8
    2Sz = 0

| **大まかなルールは次のとおりです。**

*  各行にはひと組ずつキーワード(\ ``=``\ の前)と
   パラメータ(\ ``=``\ の後)が書かれており間は\ ``=``\ で区切られています。

*  各キーワードは順不同に記述できます。

*  空白行、または\ ``//``\ で始まる行(コメントアウト)は読み飛ばされます。

*  各キーワード、パラメータの大文字\ :math:`\cdot`\ 小文字は区別されません。
   ダブルクオート、空白は無視されます。

.. *  必ず指定しなければいけないパラメータ、
   指定しない場合デフォルト値が使われるパラメータ、
   (他のパラメータの組み合わせによっては)使われないパラメータが存在します。
   使われないパラメータが指定された場合にはプログラムは終了し、
   入力ファイルをチェックするようにというメッセージが表示されます。

*  | パラメータには三種類あります。
   | 1. 必ず指定しなければいけないパラメータ (もし存在しない場合には、 :math:`{\mathcal H}\Phi` のエラーメッセージが出力され、プログラムは終了します)
   | 2. 指定しない場合デフォルト値が使われるパラメータ (もし存在しない場合は、デフォルト値が使用されます)
   | 3. 使われないパラメータ (使われないパラメータが指定された場合には"入力ファイルをチェックするように"というメッセージが表示され、プログラムは終了します。)
   | 例えば、ハイゼンベルグ模型でトランスファー積分 :math:`t` を指定した場合が相当します。 もし "model=spin" とした場合には、  "\ :math:`t`" は使用できません。

次に各キーワードの説明をします。

.. toctree::
   :maxdepth: 1

   Parameters_for_the_type_of_calculation_ja
   Parameters_for_the_lattice_ja
   Parameters_for_conserved_quantities_ja
   Parameters_for_the_Hamiltonian_ja
   Parameters_for_the_numerical_condition
   Parameters_for_the_dynamical_Greens_function        
   Parameters_for_time-evolution



