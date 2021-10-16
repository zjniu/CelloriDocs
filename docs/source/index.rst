Cellori
===================================

**Cellori** (Cell Origin) is a fast and robust intensity-based algorithm for clustered nuclei segmentation in fluorescence microscopy images. It segments nuclei by applying a Gaussian filter to smoothen out background noise, calculating local thresholds to isolate the foreground, and splitting clustered nuclei via local maxima analysis. Masks are obtained using the watershed algorithm.

Check out the :doc:`usage` section for further information, including
how to :ref:`installation` the project.

Examples
--------

<table border="1" cellpadding="5">

<tr>
<td align="center" valign="center">
Raw Image
</td>

<td align="center" valign="center">
Segmented Image
</td>
</tr>

<tr>
<td align="center" valign="center">
<img src="docs/examples/str_raw.png" width="500" alt="Raw Image (STR)" />
</td>

<td align="center" valign="center">
<img src="docs/examples/str_segmented.png" width="500" alt="Segmented Image (STR)" />
</td>
</tr>

<tr>
<td align="center" valign="center">
<img src="docs/examples/wm989_raw.png" width="500" alt="Raw Image (STR)" />
</td>

<td align="center" valign="center">
<img src="docs/examples/wm989_segmented.png" width="500" alt="Segmented Image (STR)" />
</td>
</tr>

</table>

Contents
--------

.. toctree::

   usage
   api
