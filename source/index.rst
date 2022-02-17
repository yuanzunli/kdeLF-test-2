kdeLF
=====

**kdeLF** is an MIT licensed pure-Python implementation of Yuan et al.'s
`method for estimating luminosity functions via Kernel Density Estimation (KDE) 
<https://ui.adsabs.harvard.edu/abs/2020ApJS..248....1Y>`_ and these pages will
show you how to use it.

This documentation won't teach you too much about KDE but there are a lot
of resources available for that (try `this one
<https://www.inference.org.uk/mackay/itprnn/book.html>`_).
We also `published a paper <https://arxiv.org/abs/1202.3665>`_ explaining
the kdeLF algorithm and implementation in detail.

kdeLF is being actively developed on `GitHub <https://github.com/yuanzunli/kdeLF>`_.


Basic Usage
-----------

If you want to calculate the luminosity function based on a survey data, you would do
something like:

.. code-block:: python

    import numpy as np
    import kdeLF
    from scipy import interpolate

    with open('flim.dat', 'r') as f: 
        x0, y0= np.loadtxt(f, usecols=(0,1), unpack=True)  
    f_lim = interpolate.interp1d(x0, y0,fill_value="extrapolate")
    sr=6248/((180/np.pi)**2)

    lf=kdeLF.KdeLF(sample_name='data.txt',solid_angle=sr,zbin=[0.6,0.8],f_lim=f_lim,adaptive=True)
    lf.get_optimal_h()
    lf.run_mcmc()
    lf.plot_posterior_LF(z=0.718,sigma=3)

A more complete example is available in the :ref:`quickstart` tutorial.


How to Use This Guide
---------------------

To start, you're probably going to need to follow the :ref:`install` guide to
get kdeLF installed on your computer.
After you finish that, you can probably learn most of what you need from the
tutorials listed below (you might want to start with
:ref:`quickstart` and go from there).
If you need more details about specific functionality, the User Guide below
should have what you need.

We welcome bug reports, patches, feature requests, and other comments via `the GitHub
issue tracker <https://github.com/yuanzunli/kdeLF/issues>`_, but you should check out the
`contribution guidelines <https://github.com/yuanzunli/kdeLF/blob/main/CONTRIBUTING.md>`_
first.



.. toctree::
   :maxdepth: 2
   :caption: User Guide

   user/install
   user/KdeLF
   user/faq

.. toctree::
   :maxdepth: 1
   :caption: Tutorials

   tutorials/quickstart
   tutorials/plot
   tutorials/parallel
   tutorials/KS-test
   tutorials/MCMC


License & Attribution
---------------------

Copyright 2021 Zunli Yuan and `contributors <https://github.com/yuanzunli/kdeLF/graphs/contributors>`_.

kdeLF is free software made available under the MIT License. For details
see the ``LICENSE``.

If you make use of kdeLF in your work, please cite our paper
(`arXiv <https://arxiv.org/abs/2003.13373>`_,
`ADS <https://ui.adsabs.harvard.edu/abs/2020ApJS..248....1Y>`_,
`BibTeX <https://ui.adsabs.harvard.edu/abs/2020ApJS..248....1Y/exportcitation>`_).


Changelog
---------

.. include:: ../HISTORY.rst
