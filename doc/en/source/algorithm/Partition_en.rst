.. highlight:: none

.. _Sec:sec_partion_function:

Partition function and quantities at finite-temperature
-------------------------------------------------------

To avoid overflow/underflow, we compute as follows:

Partition function

.. math::

    Z(T) &= \sum_{i=1}^N \exp\left(-\frac{E_i}{T}\right)
    \nonumber \\
    &= \exp\left(-\frac{E_1}{T}\right) \left[
    1 + \exp\left(-\frac{E_2-E_1}{T}\right)+ \exp\left(-\frac{E_3-E_1}{T}\right)
    \cdots
    + \exp\left(-\frac{E_N-E_1}{T}\right)
    \right]
    \nonumber \\
    &= \exp\left(-\frac{E_1}{T}\right) \left[
    1 + \exp\left(-\frac{E_2-E_1}{T}\right)\left[
    1 + \exp\left(-\frac{E_3-E_2}{T}\right)\left[
    1 + \dots
    \left[
    1 + \exp\left(-\frac{E_N-E_{N-1}}{T}\right)
    \right]
    \right]
    \right]
    \right]

Quantity at finite tempearture

.. math::

    O(T) &= \frac{1}{Z(T)}\sum_i O_i \exp\left(-\frac{E_i}{T}\right)
    \nonumber \\
    &= \exp\left(-\frac{E_1}{T}\right) \left[
    O_1 + O_2 \exp\left(-\frac{E_2-E_1}{T}\right) + O_3\exp\left(-\frac{E_3-E_1}{T}\right)
    \cdots
    + O_N\exp\left(-\frac{E_N-E_1}{T}\right)
    \right]
    \nonumber \\
    &= \exp\left(-\frac{E_1}{T}\right) \left[
    O_1 + \exp\left(-\frac{E_2-E_1}{T}\right)\left[
    O_2 + \exp\left(-\frac{E_3-E_2}{T}\right)\left[
    O_3 + \dots
    \left[
    O_{N-1} + O_N\exp\left(-\frac{E_N-E_{N-1}}{T}\right)
    \right]
    \right]
    \right]
    \right]

