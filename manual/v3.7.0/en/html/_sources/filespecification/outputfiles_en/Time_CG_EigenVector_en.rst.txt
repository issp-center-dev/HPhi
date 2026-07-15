.. highlight:: none

.. _Subsec:timecgeigenv:

Time_CG_EigenVector.dat
-----------------------

(For the Lanczos method) The process for calculating the eigenvector by
the CG method is outputted. An example of the file format is as follows.

::

    allocate succeed !!! 
    b[4341]=1.000000 bnorm== 1.000000 
    i_itr=0 itr=5 0.0411202543 0.0000100000 
    ...
    i_itr=0 itr=155 0.0000000058 0.0000100000 
    CG OK:   t_itr=155 
    i_itr=0 itr=155 time=0.000000  
    fabs(fabs(xb)-1.0)=0.9955114473313577
    b[4341]=0.004489 bnorm== 1.000000 
    i_itr=1 itr=5 13.0033983157 0.0000100000 
    ...
    CG OK:   t_itr=275 
    i_itr=1 itr=120 time=0.000000  
    fabs(fabs(xb)-1.0)=0.0000000000001295
    number of iterations in inv1:i_itr=1 itr=120 
    t_itr=275 0.000000

.. _file_name_3:

File name
~~~~~~~~~

*  ##_Time_CG_EigenVector.dat

## indicates a header defined by [string02] in a ModPara file.

.. raw:: latex

   \newpage