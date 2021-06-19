.. highlight:: none

.. _Subsec:timecgeigenv:

Time\_CG\_EigenVector.dat
~~~~~~~~~~~~~~~~~~~~~~~~~

(Lanczos法のみ) CG法で固有ベクトルを計算する際のログを出力します。

::

    allocate succeed !!! 
    b[4341]=1.000000 bnorm== 1.000000 
    i_itr=0 itr=5 0.0411202543 0.0000100000 
    …
    i_itr=0 itr=155 0.0000000058 0.0000100000 
    CG OK:   t_itr=155 
    i_itr=0 itr=155 time=0.000000  
    fabs(fabs(xb)-1.0)=0.9955114473313577
    b[4341]=0.004489 bnorm== 1.000000 
    i_itr=1 itr=5 13.0033983157 0.0000100000 
    …
    CG OK:   t_itr=275 
    i_itr=1 itr=120 time=0.000000  
    fabs(fabs(xb)-1.0)=0.0000000000001295
    number of iterations in inv1:i_itr=1 itr=120 
    t_itr=275 0.000000

ファイル名
^^^^^^^^^^

-  ##\_Time\_CG\_EigenVector.dat

##はModParaファイル内の[string02]で指定されるヘッダを表します。

.. raw:: latex

   \newpage