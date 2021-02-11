
# Contributions

All contributions are welcome. EPOCH is designed to be flexible and
general purpose, so changes should not affect code correctness
in the general case. For this reason most capability is in the form of
physics packages and a minimal number of conditional calls.

Wherever possible we maintain "exact" matching
of answers between minor versions, and strict backwards compatibility.

EPOCH development and maintenance
is ongoing, so please search for open issues at
https://github.com/Warwick-Plasma/epoch/issues
and consider creating an issue
with the enhancement tag before starting anything to avoid
duplication of effort. The development team is quite busy, so please allow
time for responses and code review.


# How to Contribute

In general, contributions should be created as branches and pushed to either
the main EPOCH repository or to a private fork made available to the
developers. In the former case the branch name should begin with your name.
**Branches in the main repository will be deleted after merge** so if you
wish to preserve your branch, it is best to work on a fork.
Once your branch is ready, rebase your changes onto the latest version of
of the devel branch and create a merge request.
**Submit all merge requests against the devel branch.**

If your contribution adds facility to the input deck, provide a short example
in the merge request. If it fixes issues or bugs, please reference them
by id (such as "Fixes issue #1591"). This will enable the issue to be
automatically closed once the branch has been merged.

Note that by submitting your changes for inclusion into EPOCH, you are
implicitly accepting that these will be included under the same licensing
terms as the rest of the code (GPL v3).


# Full Access

For developer access to the main EPOCH repositories, please contact
<cfsa.gitlab+admin@gmail.com>.
For small fixes patches can be offered via the issue tracker which can
be found at https://github.com/Warwick-Plasma/epoch/issues


# Coding Style

This codebase conforms to a single, consistent coding style which is
outlined in the CODING_STYLE document in the root directory.
New submissions will be rejected if they
do not conform to this style.

In particular, please note that all code must conform to the Fortran F95
standard (ISO/IEC 1539-1: 1997), and keep all lines, including comments, to
less than 80 characters long.

Please also follow the standard commit-message format.
The first line is the subject, and should generally be less than 50 characters.
The second line must be blank. Any text here is ignored.
The subsequent lines are the message body, and should generally be less
than 72 characters. You can use as many lines as you like, but be concise.


# Physics Packages

EPOCH is designed to be modular and extensible. Adding physics capability
should be done via a physics package in the physics_packages subdirectory
which can then be used
in the main code. Conditional compilation (#ifdef) should be used
**only if absolutely necessary**, such as for speed critical sections of code.
Where possible, all control should be via the input deck. Backwards
compatibility must be preserved.

Documentation will be required in Latex format.
This should cover any changes to the input deck blocks or keys, any general
limitations of the code or methods and if possible links to any methods
papers or description.
