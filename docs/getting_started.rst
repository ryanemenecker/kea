Getting Started
===============

Installation
------------

.. code-block:: bash

    pip install git+https://git@github.com/ryanemenecker/kea.git

Basic usage
-----------

The main functionality of kea is provided through the ``build_library``
function, which turns one or more protein sequences into codon-optimized
DNA sequences.

.. code-block:: python

    from kea import build_library

    # A single protein sequence
    results = build_library("MKKFLVLLFCWAVLCEHN", "yeast")

    # The returned objects expose the optimized sequences
    print(results[0].coding_sequence)
    print(results[0].gc_content_coding_sequence)

See the :doc:`user_guide` for a more detailed walkthrough.
