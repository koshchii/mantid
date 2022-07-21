.. _dns_simulation_reflection_table-ref:

Reflection Table
----------------

.. image::  ../../../images/DNS_interface_powder_tof_options_tab.png
   :align: center
   :height: 400px
\

The table of reflections can be sorted by clicking on the header labels.
"h", "k" and "l" columns contain information on the Miller indices. "F" denotes
the Lorentz corrected structure factor. It is calculated using the positions
of the atoms only if a ``.cif`` file was loaded, otherwise :math:`F^2` will
be assigned to unity. "M" denotes the multiplicity of the reflection and by
hovering over the corresponding value in the **Reflection Table** one can see
equivalent reflections. The "det_rot" column gives the lowest possible values
for detector rotation, whereas the "channel" column specifies the corresponding
detector number.

If the **Fix Omega Offset** checkbox is not selected, the values for the sample
rotation angle will be evaluated with the assumpion that **hkl1** is being parallel
to the beam if the goniometer is at 0 and dispayed in the "sample_rot" column.

To determine the omega offset, one should specify values for **det_rot**,
**det_number** and **sample_rot** of an observed reflection in the **Identify**
groupbox (**hkl1** needs to be set and **Fix Omega Offset** has to be unselected).
After that, one should press the **Calculate** button and in the "sample_rot" column
of the **Reflection Table** tab double-click on the value of the angle of a matching
reflection. This will trigger the evaluation of the omega offset and the corresponding
value will be displayed in the **Omega Offset** field. If the **Fix Omega Offset** box
is checked after that, the values in the "sample_rot" column will be updated
correspondingly with the corrected values of sample rotation angle (double-clicking
on the sample rotation angle will have no effect now).
