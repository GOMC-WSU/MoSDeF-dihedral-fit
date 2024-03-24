============
Installation
============

.. note::
    The GOMC >= v2.75 software need to be installed manually, outside of this Python install,
    with it's directory/path specified in the dihedral fit function.

Installation with `mamba <https://github.com/mamba-org/mamba>`_ (Recommended)
-----------------------------------------------------------------------------
::

    $ mamba install -c conda-forge mosdef-dihedral-fit

Install with `conda <https://repo.anaconda.com/miniconda/>`_
------------------------------------------------------------
::

    $ conda install -c conda-forge mosdef-dihedral-fit


Install an editable version via the source code
-----------------------------------------------

It is common practice to utilize a pre-packaged ``Python`` distribution like
`Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_ to
ensure all of the dependencies are installed::

    $ git clone https://github.com/GOMC-WSU/MoSDeF-dihedral-fit
    $ cd mosdef_dihedral_fit
    $ conda env create -f environment.yml
    $ conda activate mosdef_dihedral_fit
    $ pip install -e .

.. note::
    The installation instructions are for ``Unix`` and ``OSX``.
    If ``Windows`` is being used, you need to use a virtual machine or the ``Linux`` subsystem,
    since some parts of this software or its dependencies could not be compatible with ``Windows``.


Install pre-commit
------------------

To maintain uniform coding, this software utilizes the `pre-commit <https://pre-commit.com/>`_ package.

To check all the files using pre-commit, run::

     $ pre-commit run --all-files


Supported Python Versions
-------------------------

``Python 3.10 and 3.11`` are currently the only officially supported and tested version during the
software development. Older versions of ``Python`` may work, but they are not guaranteed to work.

Testing your installation
-------------------------

The ``MoSDeF-dihedral-fit`` software uses `pytest <https://docs.pytest.org/en/stable/>`_ to analyze the code for
errors, bugs, code changes, accuracy, and to verify that the existing implementation is correct.
This software auto-installs the ``pytest`` package with the ``mosdef_gomc`` environment.

To conduct these unit tests via ``pytest``, perform the following command from the base directory::

    $ pytest -v

Building the documentation
--------------------------

``MoSDeF-dihedral`` fit documentation was constructed using `sphinx <https://www.sphinx-doc.org/en/master/index.html>`_.
This software auto-installs the ``sphinx`` package with the ``mosdef_gomc`` environment.

The ``docs`` can be built locally with the following commands when in the ``docs`` directory::

    $ make html
