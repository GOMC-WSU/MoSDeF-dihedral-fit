"""General module for IO operations in MOSDEF-dihedral-fit."""

import os

from pkg_resources import resource_filename


def get_mosdef_dihedral_fit_fn(filename):
    """Get the whole path name for the file in the mosdef_dihedral_fit/utils/files directory.
    Parameters
    ----------
    filename : str
        The name of the file in the selected directory (mosdef_dihedral_fit/utils/files).
    Returns
    -------
    full_path_and_filename : str
        Full path to filename
    """
    full_path_and_filename = resource_filename(
        "mosdef_dihedral_fit", os.path.join("utils", "files", filename)
    )
    if not os.path.exists(full_path_and_filename):
        raise ValueError(
            f"ERROR: The {full_path_and_filename} does not exists."
        )
    return full_path_and_filename
