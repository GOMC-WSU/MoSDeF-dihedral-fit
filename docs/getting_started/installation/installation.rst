============
Installation
============

.. note::
    The GOMC software need to be installed manually, outside of this Python install,
    with it's directory/path specified in the dihedral fit function.

Recommended installation is with `mamba <https://github.com/mamba-org/mamba>`_
------------------------------------------------------------------------------
::

    $ mamba install -c conda-forge mosdef-dihedral-fit

Install with `conda <https://repo.anaconda.com/miniconda/>`_
------------------------------------------------------------
::

    $ conda install -c conda-forge mosdef-dihedral-fit

There is an issue building MoSDeF-GOMC version 1.0.0 with ``conda`` or ``conda-forge``
not extracting the latest ``conda`` build version. Therefore, the user can conduct
the additional command below or install using ``mamba`` because ``mamba`` is using the correct build.::

    $ conda install -c conda-forge sympy=1.10 garnett gsd pycifrw


Install an editable version from the source code
------------------------------------------------

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

To maintain uniform coding, this sofware utilizes the `pre-commit <https://pre-commit.com/>`_ package.

To check all the files, you can run::

     $ pre-commit run --all-files


Supported Python Versions
-------------------------

``Python 3.9`` is currently the only officially supported and tested version during the
software development. Older versions of ``Python`` may work, but they are not guaranteed to work.

Testing your installation
-------------------------

The ``MoSDeF-dihedral-fit`` sofware uses `pytest <https://docs.pytest.org/en/stable/>`_ to analyze the code for
errors, bugs, code changes, accuracy, and to verify that the existing implementation is correct.
This sofware auto-installs the ``pytest`` package with the ``mosdef_gomc`` environment.

To conduct these unit tests via ``pytest``, perform the following command from the base directory::

    $ pytest -v

Building the documentation
--------------------------

``MoSDeF-dihedral`` fit documentation was constructed using `sphinx <https://www.sphinx-doc.org/en/master/index.html>`_.
This sofware auto-installs the ``sphinx`` package with the ``mosdef_gomc`` environment.

The ``docs`` can be built locally with the following commands when in the ``docs`` directory::

    $ make html
