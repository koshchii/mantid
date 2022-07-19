.. _dns_elastic_sc-ref:

DNS Single Crystal elastic
==================

The single crystal elastic mode is for the reduction of polarized and unpolarized
single crystal diffraction data at the DNS instrument at MLZ.



Path / Data / Script Generator Tabs
-----------------------------------
The **path**, **data** and **script generator** tabs are described in
:ref:`DNS reduction <DNSReduction-ref>`.


Option Tab
----------
.. image::  ../images/dns_elastic_sc_options.jpg
   :height: 200px

**Wavelength** defines the used wavelength, if **get from data** is checked, it
is read from the data files. In case the selected data files have
different wavelength or the wavelength does not match the selector speed an
error is given and **get from data** is turned off. In this case, the wavelength
has to be given manually.

**Normalisation** selects whether the dataset is normalized on counting time or
monitor counts.

**Detector efficiency correction** corrects for different detector efficiency,
 angular coverage and Lorentz factor by  normalisation to vanadium data.
It is further possible to **Sum Vanadium over detector positions**, sum SF and
NSF channels and ignore the fields in the vanadium files.

**Flipping ratio Correction** corrects for the finite flipping ratio,
using NiCr data.

**Subtract instrument background from sample** allows substraction of the instrument background from the sample, the background can be scaled by a factor.

**Orientation**
is used to specify the horizontal scattering plane.
For this two vectors hklx and hkly lying in the horizontal plane must be given, they do not need to be perpendicular. In case they are not perpenciular the GUI use instead of hkly the hector which is perpedicular to hklx and lies in the plane defined by hklx and hkly.
Omegaoffset gives the deviation of the Omega rotation angle from 0 for the position when hklx is parallel to the direct beam.

**Lattice parameters**
Here the lattice parameters of the sample can be given, from which d-spacing for hklx and hkly is automatically calculated.
Which will be used to transform the thwotheta-omega space to hklx-hkly Q space.
Alternatively dx and dy, which are the d-spacing of hklx and hkly respectively, can be given manually if the corresponding checkbox is checked.

The **Binning** is normally automatically chosen to match the physical binning of the data defined by the detector bank positions an omega scan range, but can be overwritten for special cases.


Plotting Tab
------------
The **Plotting** tab offers basic plotting functionality for reduced data sets.

.. image::  ../images/dns_elastic_sc_plot.jpg
   :height: 200px

On the left a list with the different channels of a reduced datasets is shown,
which can be checked for plotting.
With the up and down arrows, one can move through single plots.

The buttons to the right have the following functionality.
1. toggle through different grid styles
2. use crystalograhic axes, this will switch between a) orthogonal axes defined by hklx and hkly or b) the crystalographic principal axes lying in the horizontal scattering plane, the latter only works, if only two principle axes change in the horizontal plane.
3. toggle on and off projections of the intensity averaged along hklx and hkly
4. toggle through drawing of borders of the quadliterals or triangles of the plot. a) no border, b) white border c) black border (might have to zoom in to see not only borders) In **Scatter** mode the button can be used to change the pointsize.
5. export the shown plot as an ascii file to the directory defined for exports in the Path tab.

Matplotlib controls:
Home button: unzooms the plot
rest as in matplotlib

On the right of the matplotlib controls, the position of the mouse cursor on the plot is shown, both in x-y coordinates, as well as in hkl. Also the Intensity with errorbar of the closest measured point is given. (This does not give the intensity of the quadliteral which could involve interpolation.)

The fields X, Y and Z below this button, allow to manually specific a region to zoom into.
The syntax is like python range (or dnsplot), so **0:2** would be from 0 to 2, $0:$ would be from 0 to the maximum, pressing enter inside will update the plot.

**log** toggles on and off logarithmic intensity scale
The dropdown list allows to select different colormaps. The button behind the dropdown list inverts the colormap color scale.
And the Fontsize dialog allows to change the fonsize, **Enter** has to be pressed inside for the change to have an effect.

The **Options** menu allows to **Change Omegaoffset**, this is a finder dialog where the user can interactively change the omegaoffset of the actual plot to search for the correct value. Pressing the arrows on the inputbox will change by plus or minus 10°.
The changes are only temporal and not applied to all plot, the user should, after the correct Offset is found, change it in the **options** tab and rerun the script.

**Change dx/dy* allows to adjust the d-spacing of hklx and hkly, again the changes are only temporal.

The **View** Menu offers control over the plot style.
**Plot Type** can be changed between triangulation where, the map is created out of triangles build between the measured points by Delaunay-triangulation. This mode works with any kind of data distribution, but is the slowest.
**Quadmesh** builds the map out of quadliterals, this is much faster and is used by the old software dnsplot, for this the number of data points in omega for every two theta position must be the same, otherwise quadliterals are not well defined.
**Scatter** gives just a plot of data points.

**Axes** allows to change the plot axes between hklx and hkly as given in the options or qx and qx, or twotheta and omega.
**Fix aspect ratio** will fix the aspect ration between the x and y axis of the plot in Q, this is especially usefull if crystallographic axes are used since then the shown angles will be correct.

**Switch axes** will switch x and y axes in the plot.

**Interpolation** allows to add additional quadliterals/triangles for a smoother looking plot. 1->9, replaces for example each qualiteral with 9 quadliterals and so on. This is a simular degree of interpolation as used by dnsplot.

**Gourad Shading**, can be toggled on and off, if it is off each quadliteral will have a solid color, determined by the average of the intensity of the datapoints on its corners. If turned on, the color will be smeared out inside the qualiteral towards its corners, giving a smoother change of color.

**Synchronize zooming** allows to synchronize the zoom betweeen different plots, if turned on for specific axis, zooming in with the mouse will be synchronous through all plots in the list.





Used By
^^^^^^^

:ref:`DNS reduction <DNSReduction-ref>`
