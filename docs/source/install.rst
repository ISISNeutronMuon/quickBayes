.. _install:

Getting Started
===============



Installation
------------

The quickBayes package is easy to install using :code:`pip`.
The library is available on `PyPi <https://pypi.org/project/quickBayes/#description/>`_ and can be installed with pip:

.. code-block:: python

    python -m pip install quickBayes

The default install does not include `gofit <https://ralna.github.io/GOFit/_build/html/index.html>`_.
To install quickBayes with gofit the command is

.. code-block:: python

    python -m pip install quickBayes[gofit]


Running the Tests
-----------------

The tests are split into three sections:

- shared
- default
- gofit

and these distinctions exist to allow the testing for both with and without gofit.
The tests in shared correspond to those that should be the same regardless of the installed packages.
The default tests are what should happen if an optional dependency is missing.
Typically this is producing an error at initialisation.
The gofit tests are to test the functionality of gofit.

To run the default tests

.. code-block:: python

    pytest test/default

and to run both the gofit and shared tests

.. code-block:: python

    pytest test/shared test/gofit


Reporting Issues
----------------

To report an issue please create a `github issue <https://github.com/ISISNeutronMuon/quickBayes/issues/>`_ with the details of the bug.


Support
-------

Support will be provided via `github issues <https://github.com/ISISNeutronMuon/quickBayes/issues/>`_, please outline the problem in a new issue.
