import os
import subprocess

from warnings import warn

def delete_directory(directory_path_and_name):
    """Delete an existing directory.

    This function deletes the directory that it was provided via the
    directory path and name

    Parameters
    ----------
    directory_path_and_name: str
        The directories' full or relative path and name of the directory to delete.

    Returns
    -------
    Deletes the selected directory, if it exists, via the terminal command
    'rm -r directory_path_and_name'
    """
    if os.path.isdir(directory_path_and_name) is True:
        delete_directory_command = f"rm -r {directory_path_and_name}"

        exec_delete_directory_command = subprocess.Popen(
            delete_directory_command, shell=True, stderr=subprocess.STDOUT
        )

        os.waitpid(
            exec_delete_directory_command.pid, os.WSTOPPED
        )  # pauses python until exec_delete_directory_command is done


def create_directory(directory_path_and_name):
    """Creates a directory if it does not currently exist

    This function creates a directory that it was provided via the
    directory path and name, if it does not currently exist

    Parameters
    ----------
    directory_path_and_name: str
        The directories full or relative path and name to create.

    Returns
    -------
    Creates the selected directory, if it does not currently exist,
    via the terminal command 'mkdir directory_path_and_name'
    """
    if os.path.isdir(directory_path_and_name) is False:
        create_directory_command = f"mkdir {directory_path_and_name}"

        exec_create_directory_command = subprocess.Popen(
            create_directory_command, shell=True, stderr=subprocess.STDOUT
        )

        os.waitpid(
            exec_create_directory_command.pid, os.WSTOPPED
        )  # pauses python until exec_create_directory_command is done

    else:
        raise ValueError(f"ERROR: The {directory_path_and_name} is trying to be created"
                         "with the 'create_directory' function, but it can not because the"
                         "it already exists.")

def delete_file(file_path_and_name):
    """Delete an existing file.

    This function deletes the file that it was provided via the
    file_path_and_name

    Parameters
    ----------
    file_path_and_name: str
        The files' full or relative path and name of the file to delete.

    Returns
    -------
    Deletes the selected file, if it exists, via the terminal command
    'rm file_path_and_name'
    """
    if os.path.isfile(file_path_and_name) is True:
        delete_file_command = f"rm {file_path_and_name}"

        exec_delete_file_command = subprocess.Popen(
            delete_file_command, shell=True, stderr=subprocess.STDOUT
        )

        os.waitpid(
            exec_delete_file_command.pid, os.WSTOPPED
        )  # pauses python until exec_delete_file_command is done