"""General module for IO operations in MOSDEF-dihedral-fit."""
import os

from pkg_resources import resource_filename


def get_mosdef_dihedral_fit_fn(filename):
    """Get the full path to one of the reference testing files provided with utils.
    In this source distribution, these files are in ``mosdef_dihedral_fit/utils/files``,
    but on installation, they're moved to somewhere in the user's python
    site-packages directory.
    Parameters
    ----------
    filename : str
        Name of the file to load (with respect to the files/folder).
    Returns
    -------
    fn : str
        Full path to filename
    """
    fn = resource_filename(
        "mosdef_dihedral_fit", os.path.join("utils", "files", filename)
    )
    if not os.path.exists(fn):
        raise IOError("Sorry! {} does not exists.".format(fn))
    return fn
