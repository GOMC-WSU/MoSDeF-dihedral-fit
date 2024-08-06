================================
How to Contribute to the Project
=================================

We welcome contributions to our project, whether it's adding a feature, fixing a bug, or improving documentation. The process is straightforward:

1. **Open an Issue on GitHub**: Start by opening an issue on GitHub to describe the change you want to make. This allows for discussion and ensures everyone is aligned before significant work is undertaken. For bug fixes, we will confirm the behavior of the bug and validate the proposed fix. For new features, we will assess the feature's relevance and discuss possible designs.  If you can provide code a short code section to show how the bug occurs or an example of the feature that you want, it will assist us in finding the bug or knowing what should be added.  For added security measures, we prefer that you provide the short code sections in the text of the GitHub Issue and not by uploaded files or compressed folders/directories.

2. **Create a Pull Request (PR)**: Once there is consensus, [create a pull request](https://help.github.com/en/articles/about-pull-requests) with your code changes. For larger features, you can create a pull request before completing the implementation to receive early feedback. Indicate that it is a work in progress by including "WIP" at the start of the PR title.

3. **Documentation and Unit Tests**: For new features, please provide sufficient documentation and unit tests:
   - **Documentation**: Include a summary of the method, explanations of all input parameters, and the expected output. A minimal example can be found [here](https://github.com/GOMC-WSU/MoSDeF-dihedral-fit/blob/main/mosdef_dihedral_fit/utils/io.py).
   - **Unit Tests**: Unit tests for the all possible options that a user can select are mandatory, as this ensures that all the code is covered and checked regularly for errors.  We use [pytest](https://docs.pytest.org/en/latest/), which can be run by executing `pytest` from the root directory of the package. These tests are automatically run as part of our CI workflow each time a change is added.

4. **Review and Merge**: Once the PR is marked as ready and all checks have passed, core developers/maintainers will review it and may suggest changes. After everyone is satisfied with the code, the pull request will be merged.

Congratulations on a successful PR, and thank you for your contribution!
