Usage
=====

.. _installation:

Installation
------------

Install Cellori from `PyPI <https://pypi.org/project/cellori/>`_.

.. code-block::

    pip install cellori
    
.. _tutorial:

Guided Segmentation Tutorial
----------------

1. Download the `sample image <https://github.com/SydShafferLab/Cellori/tree/main/docs/demo/wm989.tif>`_ (WM989 cells).
2. Create a Python script with the following lines of code. Make sure that wm989.tif is in your current working directory.

.. code-block:: python

    from cellori import Cellori

    Cellori('wm989.tif').gui()
    
3. After running the script, you will be greeted with the following window. Notice that no cells are currently in view.

.. image:: ../demo/gui1.png
           :width: 1000
           :alt: GUI 1
           
4. Change the preview region (indicated by the red box) by clicking anywhere on the left panel or using your arrow keys. The center of each nucleus is marked with a red dot, and the total count in the preview region is shown above the right panel. Here, we are looking at a colony of cells in the bottom left of the image.

.. image:: ../demo/gui.png
           :width: 1000
           :alt: GUI 2
           
5. Automatic parameter detection should have already chosen values that work well, but they can be manually adjusted if desired. Here is a brief description of each parameter.

    * Sigma: Gaussian sigma used for background denoising.
    * Block Size: Odd size of pixel neighborhood which is used for local thresholding (e.g., 3, 5, 7, ..., 21).
    * Nuclei Diameter: Estimated lower bound of nuclei diameters. Any objects smaller than this threshold will not be considered for segmentation.

6. We will first explore the efforts of the sigma parameter. A higher sigma results in more blurring, which reduces the issues of background noise and over-segmentation of single nuclei. If we look at segmentation when sigma is 0.5, notice that some single nuclei are being split up into two or even three separate nuclei. However, a sigma that is too high could fail to split clustered nuclei, or worse, miss nuclei altogether, as seen in the segmentation when sigma is 3.5.

.. list-table::
   :widths: 33 33 33
   :header-rows: 1

   * - Sigma = 0.5
     - Sigma = 1.5
     - Sigma = 3.5
   * - .. image:: ../demo/sigma0.5.png
           :width: 300
           :alt: Sigma = 0.5
     - .. image:: ../demo/default.png
           :width: 300
           :alt: Sigma = 1.5
     - .. image:: ../demo/sigma3.5.png
           :width: 300
           :alt: Sigma = 3.5

7. Click on the "Segment" button to segment the entire image.
