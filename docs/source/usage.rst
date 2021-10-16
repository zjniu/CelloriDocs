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
.. image:: ../examples/gui1.png
           :width: 500
           :alt: GUI 1
4. Change the preview region (indicated by the red box) by clicking anywhere on the left panel or using your arrow keys. The center of each nucleus is marked with a red dot, and the total count in the preview region is shown above the right panel.

