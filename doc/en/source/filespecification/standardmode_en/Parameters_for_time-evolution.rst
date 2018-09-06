.. highlight:: none

Parameters for time-evolution
-----------------------------

*  ``dt``

   **Type :** Positive Double(Default value : ``0.1``)

   **Description :** The width of time steps.

*  ``PumpType``

   **Type :** String (Chosen from ``"Quench"``, ``"Pulse Laser"``,
   ``"AC Laser"``, and ``"DC Laser"``.Default value : ``"Quench"``)

   **Description :** The type of time-dependent Hamiltonian. For
   ``"Quench"``, two body operator
   :math:`U_{\rm quench} \sum_i n_{i \uparrow} n_{i \downarrow}` is
   added. For ``"Pulse Laser"``, ``"AC Laser"``, and ``"DC Laser"``, the
   hopping term is modulated as
   :math:`-\sum_{i j \sigma} t_{i j} \exp[-{\bf A}(t) \cdot ({\bf R}_i-{\bf R}_j)/(2\pi)] c_{i \sigma} c_{j \sigma}`,
   where :math:`{\bf A}(t)` is the vector potential which is defined as
   :math:`{\bf A}(t) = {\bf A}_0 \exp[-(t-t_0)^2/(2 t_{\rm dump}^2)] \cos[\omega (t-t_0)]`,
   :math:`{\bf A}(t) = {\bf A}_0 \sin[\omega (t-t_0)]`, and
   :math:`{\bf A}(t) = {\bf A}_0 t` for for ``"Pulse Laser"``,
   ``"AC Laser"``, and ``"DC Laser"``, respectively.

   ``potential.dat`` file is written for displaying the vector potential
   and the electrical field at each time step.

*  ``Uquench``

   **Type :** Double(Default value : ``0.0``)

   **Description :** :math:`U_{\rm quench}`

*  ``freq``

   **Type :** Double(Default value : ``0.1``)

   **Description :** :math:`\omega`

*  ``tshift``

   **Type :** Double(Default value : ``0.0``)

   **Description :** :math:`t_0`

*  ``tdump``

   **Type :** Double (Default value : ``0.1``)

   **Description :** :math:`t_{\rm dump}`

*  ``VecPotW``, ``VecPotL``

   **Type :** Double (Default value : ``0.0``)

   **Description :** The vector potential (fractional coordinate of the
   reciprocal space) at :math:`t=t_0`. The reciprocal lattice vector is
   computed from the direct lattice vector shown in
   :numref:`fig_chap04_1_lattice` , :numref:`fig_chap04_1_honeycomb` ,
   :numref:`fig_ladder` , :numref:`fig_kagome` .

.. raw:: latex

   \newpage
