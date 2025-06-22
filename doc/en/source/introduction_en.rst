.. highlight:: none

*********************************
What is :math:`{\mathcal H}\Phi`?
*********************************

What is :math:`{\mathcal H}\Phi`?
=================================

Comparison between experimental observation and theoretical analysis is a crucial step in condensed-matter physics research. The temperature dependence of specific heat and magnetic susceptibility, for example, has been studied to extract the nature of low energy excitations of and magnetic interactions between electrons, respectively, through comparison with theories such as Landau's Fermi liquid theory and the Curie-Weiss law.

For the flexible and quantitative comparison of theoretical and  experimental data, the exact diagonalization approach [1]_ is one of the most reliable numerical tools that requires no approximation or inspiration of genius. For the last few decades, a numerical diagonalization package for quantum spin Hamiltonians, TITPACK, developed by Prof. Hidetoshi Nishimori of Tokyo Institute of Technology, has been widely used in the condensed-matter physics community. Nevertheless, limited computational resources have hindered the ability of non-expert users to apply the package to quantum systems with a large number of electrons or spins.

In contrast, the recent and rapid development of a parallel computing infrastructure has opened up new avenues for user-friendly larger scale diagonalizations up to 18-site Hubbard clusters or 36 :math:`S=1/2` quantum spins. In addition, recent advances in quantum statistical mechanics [2]_ [3]_ [4]_ [5]_ allow the finite temperature properties of quantum many-body systems to be calculated at computational costs similar to those of the calculations of ground state properties, which also allows theoretical results for the temperature dependence of, for example, specific heat and magnetic susceptibility, to be compared with experimental results quantitatively [6]_ . To utilize the parallel computing infrastructure with narrow bandwidth and distributed-memory architectures, efficient, user-friendly, and highly parallelized diagonalization packages are highly desirable.

:math:`{\mathcal H}\Phi`, a flexible diagonalization package for solving quantum lattice Hamiltonians, has been developed as a descendant of the pioneering package TITPACK. The Lanczos method for calculations of the ground state and a few excited states properties, as well as finite temperature calculations based on thermal pure quantum states [5]_ , are implemented in the :math:`{\mathcal H}\Phi` package, with an easy-to-use and flexible user interface. By using :math:`{\mathcal H}\Phi`, you can analyze a wide range of quantum lattice Hamiltonians including simple Hubbard and Heisenberg models, multi-band extensions of the Hubbard model, exchange couplings that break the SU(2) symmetry of quantum spins, such as Dzyaloshinskii-Moriya and Kitaev interactions, and Kondo lattice models describing itinerant electrons coupled with quantum spins. :math:`{\mathcal H}\Phi` calculates a variety of physical quantities, such as internal energy at zero temperature or finite temperatures, temperature dependence of specific heat, and charge/spin structure factors. A broad spectrum of users including experimental scientists is cordially welcome.

License
-------

The distribution of the program package and the source codes for :math:`{\mathcal H}\Phi` follow GNU General Public License version 3 (GPL v3) or later. We hope that you cite the reference, `Comp. Phys. Commun. 217 (2017) 180-192 <https://www.sciencedirect.com/science/article/pii/S0010465517301200?via%3Dihub>`_ , when you publish the results using :math:`{\mathcal H}\Phi` (hphi).

Copyright
---------

© *2015- The University of Tokyo. All rights reserved.*
This software was developed with the support of \"*Project for advancement of software usability in materials science*\" of The Institute for Solid State Physics, The University of Tokyo. 

Contributors
------------

This software was developed by the following contributors.

* ver.3.5.2 (released on 2024/3/07)

* ver.3.5.1 (released on 2022/6/14)
  
* ver.3.5 (released on 2021/9/29)

* ver.3.4 (released on 2020/6/2)

  * Developers

    * | Kota Ido
      | (The Institute for Solid State Physics, The University of Tokyo)
    * | Mitsuaki Kawamura
      | (The Institute for Solid State Physics, The University of Tokyo)
    * | Yuichi Motoyama
      | (The Institute for Solid State Physics, The University of Tokyo)
    * | Kazuyoshi Yoshimi
      | (The Institute for Solid State Physics, The University of Tokyo)
    * | Takahiro Misawa
      | (Waseda Research Institute for Science and Engineering (RISE), Waseda University)
    * | Youhei Yamaji
      | (Department of Applied Physics, The University of Tokyo)
    * | Synge Todo
      | (Department of Physics, The University of Tokyo)
    * | Yusuke Konishi
      | (Academeia Co., Ltd.)
   
  * Project coordinator

    * | Naoki Kawashima
      | (The Institute for Solid State Physics, The University of Tokyo)


* ver.3.3 (released on 2019/7/19)

* ver.3.2 (released on 2019/4/27)

* ver.3.1 (released on 2018/9/3)

* ver.3.0 (released on 2017/12/22)

  * Developers

    * | Takahiro Misawa
      | (The Institute for Solid State Physics, The University of Tokyo)
    * | Kazuyoshi Yoshimi
      | (The Institute for Solid State Physics, The University of Tokyo)
    * | Mitsuaki Kawamura
      | (The Institute for Solid State Physics, The University of Tokyo)
    * | Kota Ido
      | (Department of Applied Physics, The University of Tokyo)
    * | Youhei Yamaji
      | (Department of Applied Physics, The University of Tokyo)
    * | Synge Todo
      | (Department of Physics, The University of Tokyo)
   
  * Project coordinator

    * | Naoki Kawashima
      | (The Institute for Solid State Physics, The University of Tokyo)

* ver.2.0 (released on 2017/4/11)
* ver.1.2 (released on 2016/11/14)
* ver.1.1 (released on 2016/5/13)
* ver.1.0 (released on 2016/4/5)

  * Developers

    * | Takahiro Misawa
      | (Department of Applied Physics, The University of Tokyo)
    * | Kazuyoshi Yoshimi
      | (The Institute for Solid State Physics, The University of Tokyo)
    * | Mitsuaki Kawamura
      | (The Institute for Solid State Physics, The University of Tokyo)
    * | Youhei Yamaji
      | (Department of Applied Physics, The University of Tokyo)
    * | Synge Todo
      | (Department of Physics, The University of Tokyo)
   
  * Project coordinator

    * | Naoki Kawashima
      | (The Institute for Solid State Physics, The University of Tokyo)
   
Operating environment
=====================

:math:`{\mathcal H}\Phi` was tested on the following platforms

* The supercomputer system-B \"ohtaka\" in ISSP
* Linux PC + Intel compiler
* Linux PC + GCC.
* Mac + GCC.

.. [1] \E. Dagotto, Rev. Mod. Phys. **66**, 763-840 (1994).
.. [2] \M. Imada, M. Takahashi, Journal of the Physical Society of Japan **55**, 3354-3361 (1986).
.. [3] \J. Jaklič, P. Prelovšek, Phys. Rev. B **49**, 5065-5068 (1994).
.. [4] \A. Hams, H. De Raedt, Phys. Rev. E **62**, 4365-4377 (2000).
.. [5] \S. Sugiura, A. Shimizu, Phys. Rev. Lett. **108**, 240401 (2012).
.. [6] \Y. Yamaji, Y. Nomura, M. Kurita, R. Arita, M. Imada, Phys. Rev. Lett. **113**, 107201 (2014).
