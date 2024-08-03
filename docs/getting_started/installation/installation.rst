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

Install `GOMC <https://gomc-wsu.org/>`_
-----------------------------------------------------------------------------
::

    $ pip install cmake
    $ mkdir GOMC_INSTALL
    $ cd GOMC_INSTALL
    $ git clone https://github.com/GOMC-WSU/GOMC.git --branch v2.75a
    $ cd GOMC
    $ chmod u+x metamake.sh
    $ ./metamake.sh

.. note::
   Add GOMC your path so it is called automatically with path entereed as "". 
   To so this, you must add the following to your ".bashrc" file or 
   equivalent file (if Mac maybe ".zshrc" depending on how you set up your 
   specific system). If a different than ".bashrc", replace ".bashrc" with ".XXXrc" file, 
   where XXX is your scripting language of choice, and can be found by running "$ echo $0". The first 2 commands add a 
   space and title to the command so it is known and traceable later. 
   The commands are listed below: ::

   $ echo "" >> ~/.bashrc
   $ echo "# Add the GOMC binary files to the path" >> ~/.bashrc
   $ cwd=$(pwd) # you should still be in the GOMC directory cloned from github.
   $ echo "export PATH=${cwd}/bin${PATH:+:${PATH}}" >> ~/.bashrc

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
This software auto-installs the ``pytest`` package with the ``MoSDeF-dihedral-fit`` environment.

To conduct these unit tests via ``pytest``, perform the following command from the base directory::

    $ pytest -v

.. note::
    In the ``MoSDeF-dihedral-fit/mosdef_dihedral_fit/tests/test_fit_dihedral_with_gomc.py`` file, 
    you will have to set the ``gomc_binary_directory`` variable (near the top of the file) 
    equal to your full and local GOMC binary directory in order for the 
    ``test_fit_dihedral_with_gomc.py`` file not to fail all the unit tests (``pytest``).
    
    Example: **gomc_binary_directory = /Users/brad/Programs/GOMC/GOMC_2_75a/bin**  
    

Building the documentation
--------------------------

``MoSDeF-dihedral-fit`` documentation was constructed using `sphinx <https://www.sphinx-doc.org/en/master/index.html>`_.
The ``sphinx`` software may need to be installed separately to avoid dependency conflicts. 
If ``sphinx`` is not automatically provided, the correct ``sphinx`` package can be build after creating 
a new conda environment using the ``environment_docs.yml`` file in the ``MoSDeF-dihedral-fit/docs`` 
directory, located on ``MoSDeF-dihedral-fit`` GitHub's main repository or GitHub's releases for a specific version.

Once the correct ``sphinx`` package is installed, 
the ``docs`` can be built locally with the following commands when in the ``docs`` directory::

    $ make html
