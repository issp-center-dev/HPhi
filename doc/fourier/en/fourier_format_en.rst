.. _fileformat:

File format
===========

.. _geometry:

Geometry
--------

The file name in the :ref:`tutorial` is ``geometry.dat``.
When we use Standard mode of mVMC/:math:`{\mathcal H}\Phi`,
this file is generated automatically.
Therefore we do not have to care it.

::

   1.000000     0.000000     0.000000  (1)
   0.000000     1.000000     0.000000  (1)
   0.000000     0.000000     1.000000  (1)
   0.000000     0.000000     0.000000  (2)
   4 0 0                               (3)
   0 4 0                               (3)
   0 0 1                               (3)
   0.000000     0.000000     0.000000  (4)
   1.000000     0.000000     0.000000  (4)
   2.000000     0.000000     0.000000  (4)
   3.000000     0.000000     0.000000  (4)
   0.000000     1.000000     0.000000  (4)
   1.000000     1.000000     0.000000  (4)
   2.000000     1.000000     0.000000  (4)
   3.000000     1.000000     0.000000  (4)
   0.000000     2.000000     0.000000  (4)
   1.000000     2.000000     0.000000  (4)
   2.000000     2.000000     0.000000  (4)
   3.000000     2.000000     0.000000  (4)
   0.000000     3.000000     0.000000  (4)
   1.000000     3.000000     0.000000  (4)
   2.000000     3.000000     0.000000  (4)
   3.000000     3.000000     0.000000  (4)

#. The unit lattice vectors. Arbitrary unit.
#. The phase for the one-body term across boundaries of the simulation cell (degree unit).
#. Three integer vector specifying the shape of the simulation cell.
   They are the same as the input parameters ``a0W``, ``a0L``, ``a0H``, ``a1W``...
   in Standard mode.
#. The position of each site. The fractional coordinate is used.
   

One- and Two-body correlation function in the site representation
-----------------------------------------------------------------

.. _greenindex:

Specify the index of correlation function to be computed
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Specify the index of correlation functions
computed with mVMC/:math:`{\mathcal H}\Phi`.
When we use the standard mode, this file is generated automatically.
The general description is written in the manuals for mVMC/:math:`{\mathcal H}\Phi`.
The file names in the :ref:`tutorial` are ``greenone.def`` (one body) and ``greentwo.def`` (two body).

For calculating correlation functions in :ref:`supported`,
indices must be specified as follows:

- :math:`\langle {\hat c}_{{\bf k} \uparrow}^{\dagger} {\hat c}_{{\bf k} \uparrow}\rangle`

  :math:`\langle {\hat c}_{i \uparrow}^{\dagger} {\hat c}_{j \uparrow}\rangle`
  with :math:`(i, j)` ranging on the whole site.
  
- :math:`\langle {\hat c}_{{\bf k} \downarrow}^{\dagger} {\hat c}_{{\bf k} \downarrow}\rangle`

  :math:`\langle {\hat c}_{i \downarrow}^{\dagger} {\hat c}_{j \downarrow}\rangle`
  with :math:`(i, j)` ranging on the whole site.
  
- :math:`\langle {\hat \rho}_{\bf k} {\hat \rho}_{\bf k}\rangle` and
  :math:`\langle {\hat S}_{\bf k}^{z} {\hat S}_{\bf k}^{z} \rangle`

  :math:`\langle {\hat c}_{i \sigma}^{\dagger} {\hat c}_{i \sigma} {\hat c}_{j \sigma'}^{\dagger} {\hat c}_{j \sigma'}\rangle`
  with :math:`(i, j)` ranging on the whole site and
  :math:`(\sigma, \sigma')` ranging from :math:`\uparrow` to :math:`\downarrow`.

- :math:`\langle {\hat S}_{\bf k}^{+} {\hat S}_{\bf k}^{-} \rangle` and
  :math:`\langle {\hat {\bf S}}_{\bf k} \cdot {\hat {\bf S}}_{\bf k} \rangle`

  For :math:`{\mathcal H}\Phi`,
  :math:`\langle {\hat c}_{i \sigma}^{\dagger} {\hat c}_{i -\sigma} {\hat c}_{j -\sigma}^{\dagger} {\hat c}_{j \sigma}\rangle`
  with :math:`(i, j)` ranging on the whole site and
  :math:`\sigma` ranging from :math:`\uparrow` to :math:`\downarrow`.
  For mVMC,
  :math:`\langle {\hat c}_{i \sigma}^{\dagger} {\hat c}_{j \sigma} {\hat c}_{j -\sigma}^{\dagger} {\hat c}_{i -\sigma}\rangle`
  with :math:`(i, j)` ranging on the whole site and
  :math:`\sigma` ranging from :math:`\uparrow` to :math:`\downarrow`.
  In the both cases, please care the order of operators.
  
In the default settings of Standard mode (``outputmode="corr"``),
the above indices are specified automatically.
Therefore we do not have to care it.

.. _zvocisajs:

Results of correlation function in the site representation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The correlation functions having the indices specified in :ref:`greenindex`
are computed by mVMC/:math:`{\mathcal H}\Phi`,
and written to files.
The general description of this file is written in the manuals of mVMC/:math:`{\mathcal H}\Phi`.
File names in the :ref:`tutorial` are
``output/zvo_cisajs_001.dat`` and ``output/zvo_cisajscktalt_001.dat`` (mVMC), or
``output/zvo_cisajs.dat`` and ``output/zvo_cisajscktalt.dat`` (:math:`{\mathcal H}\Phi`).

The utility ``fourier`` reads these files before the calculation.
If some of the correlation functions with indices written in :ref:`greenindex` are lacking
(for example, because Standard mode was not used),
this utility assume them as 0.

.. _zvocorr:

Correlation functions in the primitive Brillouin zone
-----------------------------------------------------

This file contains the Fourier-transformed correlation function and
generated by the utility ``fourier``.
The file name in the :ref:`tutorial` is ``output/zvo_corr.dat``.

::
   
   #HPhi          16                                              (1)
   # kx[1] ky[2] kz[3](Cart.) UpUp[4,5] (Re. Im.) DownDown[6,7]   (2)
   # Density[8,9] SzSz[10,11] S+S-[12,13] S-S+[14,15]             (2)
   #k-offset      0.0000000      0.0000000      0.0000000         (3)
   0.00000E+00    0.00000E+00    0.00000E+00    0.31250E-01  .... (4)
   0.15708E+01    0.00000E+00    0.00000E+00    0.31250E-01  .... (4)
   :                                                               :

#. ``"#HPhi"`` for the output of ``HPhi``,
   ``"#mVMC"`` for the output of ``vmc.out``
   The subsequent integer indicate the number of :math:`k` points in the primitive Brillouine zone.
#. The description of the quantities in each column.
#. The :math:`k` offset for the one-body correlation function.
   That is to say, the one-body correlation function in the 4-7 columns are those
   at the :math:`k` point shifted from that point in the 1-3 column.
#. The :math:`k` point (Cartesian) and correlation functions.
   The real- and the imaginary-part of each correlation function are written.
   
.. _kpoint:

*k*\-point file for corplot
---------------------------

This file is generated by ``fourier`` and
read by ``corplot`` when the correlation function is plotted.
The file name is ``kpoint.dat``.

::
   
   81           9                                      (1)
   0.62832E+01    0.00000E+00    0.00000E+00           (2)
   0.00000E+00    0.62832E+01    0.00000E+00           (2)
   0.00000E+00    0.00000E+00    0.62832E+01           (2)
   -0.62832E+01   -0.62832E+01    0.00000E+00      1   (3)
   -0.47124E+01   -0.62832E+01    0.00000E+00      2   (3)
   -0.31416E+01   -0.62832E+01    0.00000E+00      3
   -0.15708E+01   -0.62832E+01    0.00000E+00      4
   0.00000E+00   -0.62832E+01    0.00000E+00      1
   0.15708E+01   -0.62832E+01    0.00000E+00      2
   0.31416E+01   -0.62832E+01    0.00000E+00      3
   0.47124E+01   -0.62832E+01    0.00000E+00      4

#. The total number of :math:`k` points plotted by ``corplot`` and
   the number of columns for displaying by splot of gnuplot.
#. Reciprocal lattice vectors (Cartesian coordinate).
#. The :math:`k` vector (Cartesian) and
   the index of the equivalent :math:`k` point in the primitive Brillouin zone.
   This number is the same as that in :ref:`zvocorr`
   
.. _gnuplot:

gnuplot script
--------------

This file is generated by ``corplot``,
and read from gnuplot launched automatically.
We also can launch gnuplot independently and ``load`` this script.
The file name is ``correlation.gp``.

.. code-block:: gnuplot

   #set terminal pdf color enhanced \    (1)
   #dashed dl 1.0 size 20.0cm, 20.0cm    (1)
   #set output 'correlation.pdf'         (1)
   #set view 60.0, 30.0                  (1)

   set view equal xy
   set ticslevel 0
   set hidden3d
   set xlabel 'kx'
   set ylabel 'ky'
   set zrange [    0.25000E-10:    0.18435E+00]

   set pm3d
   set pm3d interpolate 5, 5
   set view 0.0, 0.0

   #####  Set Brillouin-Zone Boundary  #####

   set arrow from    -0.31416E+01,   -0.31416E+01,    ...
   set arrow from    -0.31416E+01,    0.31416E+01,    ...
   :
   #####  End Set Brillouin-Zone Boundary  #####

   splot \
   'correlation.dat' u 1:2:3 w l tit '1' (2)
   pause -1

#. When we want to write the figure to a file,
   this line is uncommented.
   For pasting this figure on the paper etc.,
   we write the setting of font, line-color, and so on.
   For more details, please see the manual of gnuplot.
#. Plotting the file in :ref:`correlation`.

.. _correlation:

Correlation function at wide range of *k*
-----------------------------------------

This file is generated by ``corplot``, and
read from gnuplot through :ref:`gnuplot`.
The file name is ``correlation.dat``.

::

   -0.62832E+01   -0.62832E+01    0.18435E+00    0.00000E+00
   -0.47124E+01   -0.62832E+01    0.36159E-01    0.00000E+00
   -0.31416E+01   -0.62832E+01    0.20921E-01    0.00000E+00
   -0.15708E+01   -0.62832E+01    0.36159E-01    0.00000E+00
    0.00000E+00   -0.62832E+01    0.18435E+00    0.00000E+00
    0.15708E+01   -0.62832E+01    0.36159E-01    0.00000E+00
    0.31416E+01   -0.62832E+01    0.20921E-01    0.00000E+00
    0.47124E+01   -0.62832E+01    0.36159E-01    0.00000E+00
    0.62832E+01   -0.62832E+01    0.18435E+00    0.00000E+00

   -0.62832E+01   -0.47124E+01    0.36159E-01    0.00000E+00
   -0.47124E+01   -0.47124E+01    0.20921E-01    0.00000E+00
   -0.31416E+01   -0.47124E+01    0.11372E-01    0.00000E+00
   :

The 1st and the 2nd column contains the :math:`k` vector (Cartesian).
3rd and the 4th column contains the correlation function and its standard error, respectively.
