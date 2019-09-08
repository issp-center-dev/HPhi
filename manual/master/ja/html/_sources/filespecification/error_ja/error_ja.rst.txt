エラーメッセージ一覧
--------------------

-  ``ERROR ! Unsupported Keyword !``

   存在しないキーワードを指定した場合に表示され、プログラムは停止します。

-  ``"ERROR !  Keyword`` ``is duplicated !``

   同じキーワードを2回指定した場合に表示され、プログラムは停止します。

-  ``ERROR ! Unsupported Solver :``  *solver*

-  ``ERROR ! Unsupported Model :``  *model*

-  | ``Sorry, this system is unsupported in the STANDARD MODE...``
   | ``Please use the EXPART MODE, or write a NEW FUNCTION and post it us.``

   ``method``, ``model``, ``lattice``\ のどれかまたは複数に
   サポートしていないパラメーターを入れた場合、プログラムは停止します。

-  ``ERROR ! abs(2 * Sz) > nsite in Hubbard model !``

-  ``ERROR ! Nelec > 2 * nsite in Hubbard model !``

-  ``ERROR ! (nelec + 2 * Sz) % 2 != 0 in Hubbard model !``

-  ``ERROR ! nelec <= nsite && 2 * |Sz| > nelec in Hubbard model !``

-  ``ERROR ! nelec > nsite && 2 * |Sz| > 2 * nsite - nelec in Hubbard model !``

-  ``ERROR ! abs(2 * Sz) > nsite in Spin model !``

-  ``ERROR ! (nsite + 2 * Sz) % 2 != 0 in Spin model !``

-  ``ERROR ! abs(2 * Sz) > nsite in Hubbard model !``

-  ``ERROR ! Nelec_cond / 2 + Nelec_loc > nsite in Kondo model !``

-  ``ERROR ! (nelec_cond + nelec_loc + 2 * Sz) % 2 != 0 in Kondo model !``

-  ``ERROR ! nelec_cond <= nsite / 2 && 2 * |Sz| > nelec_cond + nelec_loc ...``

-  ``ERROR ! nelec_cond > nsite / 2 && abs(Sz2) > nsite / 2 * 3 - nelec...``

   カノニカル集団の計算において、
   入力されたサイト数、電子数、全スピン\ :math:`z`\ 成分が実現できない組み合わせである場合
   (例えば、電子数がサイト数の2倍よりも大きい、など)プログラムは停止します。

-  | ``Check ! `` `` is SPECIFIED but will NOT be USED.``
   | ``        Please COMMENT-OUT this line``
   | ``        or check this input is REALLY APPROPRIATE for your purpose !``

   使われないパラメーターを指定した時には、ユーザーに入力ファイルの確認を促しプログラムは停止します。
   実際に必要のないパラメーターの場合は該当する行を削除もしくはコメントアウトしてください。

-  ``ERROR ! ``\ `` is NOT specified !``

   必ず指定しなければならないキーワードが指定されていない場合にはプログラムは停止します。

-  ``=`` ``######  DEFAULT VALUE IS USED  ######``

   これはエラーメッセージではありません。入力ファイルで指定がなかったためにデフォルト値が使われたことを知らせるメッセージです。
