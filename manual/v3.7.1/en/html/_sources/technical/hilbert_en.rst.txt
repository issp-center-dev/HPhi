Hilbert Space Construction
==========================

This section describes how :math:`{\mathcal H}\Phi` constructs and manages the Hilbert space.

Bit Representation
------------------

:math:`{\mathcal H}\Phi` uses bit representation to efficiently store and manipulate quantum states.

**Hubbard model (4 states per site):**
  Each site has 4 possible states: empty, spin-up, spin-down, and doubly occupied.
  These are encoded using 2 bits per site:

  - ``00``: empty
  - ``01``: spin-up electron
  - ``10``: spin-down electron
  - ``11``: doubly occupied

**Spin-1/2 model (2 states per site):**
  Each site has 2 possible states: spin-up and spin-down.
  These are encoded using 1 bit per site:

  - ``0``: spin-down
  - ``1``: spin-up

**SpinlessFermion model (2 states per site):**
  Each site is either empty or occupied by one fermion:

  - ``0``: empty
  - ``1``: occupied

State Indexing
--------------

States are indexed by converting the bit representation to an integer.
For a system with :math:`N` sites, the total number of states is:

- Hubbard: :math:`4^N`
- Spin-1/2: :math:`2^N`
- SpinlessFermion: :math:`2^N`

Symmetry Restrictions
---------------------

To reduce the computational cost, :math:`{\mathcal H}\Phi` can restrict the Hilbert space
using conserved quantum numbers.

**Particle number conservation:**
  For canonical ensembles, the total number of electrons is fixed.
  Only states with the specified particle number are included.

**Sz conservation:**
  For spin systems, the total :math:`S_z` can be fixed.
  This reduces the Hilbert space dimension significantly.

**Example:**
  For a half-filled 8-site Hubbard model with :math:`S_z = 0`:

  - Full Hilbert space: :math:`4^8 = 65,536` states
  - With particle conservation (4 electrons): reduced
  - With :math:`S_z = 0` (2 up, 2 down): :math:`\binom{8}{2}^2 = 784` states

List Arrays
-----------

:math:`{\mathcal H}\Phi` uses list arrays to map between state indices and bit representations:

``list_1[i]``
  Maps state index :math:`i` to its bit representation

``list_2_1[b]``, ``list_2_2[b]``
  Inverse mapping from bit representation to state index
  (split for efficiency in large systems)

These arrays are essential for applying operators and performing MPI communication.
