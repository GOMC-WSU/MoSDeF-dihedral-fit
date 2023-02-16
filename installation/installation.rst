============
Installation
============
.. image:: https://img.shields.io/badge/license-MIT-blue.svg
    :target: http://opensource.org/licenses/MIT

Install with `mamba <https://github.com/mamba-org/mamba>`_
----------------------------------------------------------

Installing with ``mamba`` is the recommended installation method.

::

    $ mamba install -c conda-forge mosdef-dihedral-fit

Install with `conda <https://repo.anaconda.com/miniconda/>`_
------------------------------------------------------------
::

    $ conda install -c conda-forge mosdef-dihedral-fit

MoSDeF-GOMC version 1.0.0 is a dependency of this software, and there is an issue
pulling the latest ``conda`` or ``conda-forge`` build version of MoSDeF-GOMC version 1.0.0;
Therefore, the user can run the additional command below to fix this issue.
Alternatively, the user can install this software using ``mamba``, as ``mamba`` correctly pulls the
correct build version of MoSDeF-GOMC.::

    $ conda install -c conda-forge sympy=1.10 garnett gsd pycifrw


Install an editable version from the source code
------------------------------------------------

The best practice is to use a pre-packaged Python distribution like
`Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_,
which typically guarantees all of the dependencies are installed::

    $ git clone https://github.com/GOMC-WSU/MoSDeF-GOMC
    $ cd mosdef_dihedral_fit
    $ conda env create -f environment.yml
    $ conda activate mosdef_dihedral_fit
    $ pip install -e .

.. note::
    The installation instructions above are for Unix and OSX.
    If this software is being installed on Windows, a virtual machine or the Linux subsystem,
    should be used because some parts of this software and its dependencies may not be compatible with Windows.


Install pre-commit
------------------

The `pre-commit <https://pre-commit.com/>`_ packages are auto-installed with this code so it can preserve
uniform code formatting in the mosdef_dihedral_fit environment.

To use **pre-commit** to check all the files, run the following command::

     $ pre-commit run --all-files


Supported Python Versions
-------------------------

Python 3.9 is the officially supported and tested version during the development of the distributed final product.
Other versions of Python may work, but the developers can not guarantee they will work.

Testing your installation
-------------------------

MoSDeF-dihedral-fit periodically tests the code for accuracy, bug, code changes and updates, and also tests
the existing code is correct via `pytest <https://docs.pytest.org/en/stable/>`_.
The mosdef_dihedral_fit environment automatically installs the pytest package.

The user can manually run these unit tests from the base directory using the following command::

    $ pytest -v

Building the documentation
--------------------------

MoSDeF-dihedral-fit documentation is constructed via the
`sphinx <https://www.sphinx-doc.org/en/master/index.html>`_ software.
The mosdef_dihedral_fit environment automatically installs the sphinx package.

The documentation (``docs``) can be built locally from the ``docs`` directory using the following commands.::

    $ make html
