User Guide
===============

This page details how to use kea.

Building a library
------------------

``build_library`` accepts a single protein sequence (``str``), a list of
sequences, or a ``{name: sequence}`` dictionary, together with a codon
frequency table (one of ``"yeast"``, ``"s288c"``, ``"s288c_unweighted"`` or a
custom ``{amino_acid: {codon: frequency}}`` dict).

.. code-block:: python

    from kea import build_library

    sequences = {
        "protein1": "MKKFLVLLFCWAVLCEHN",
        "protein2": "MVLSEGEWQLVLHVWAKV",
    }

    results = build_library(
        sequences,
        "yeast",
        target_gc_range=(0.4, 0.6),
        total_length=120,
        adapter_5_prime="GGTCTC",
        adapter_3_prime="GAGACC",
    )

Each item in ``results`` is a ``Sequence`` object exposing, among others,
``protein_sequence``, ``coding_sequence``, ``full_dna_sequence``,
``gc_content_coding_sequence``, ``gc_content_full_sequence``,
``correct_coding_translation`` and ``correct_full_translation``.

Controlling GC content
----------------------

When ``target_gc_range`` is provided, kea fine-tunes the optimized sequence so
its GC content falls within the requested range (within ``gc_tolerance``). If
no range is provided, GC content is not optimized.

Padding to a fixed length
-------------------------

Set ``total_length`` to pad every sequence to a common length. The
``pad_location`` argument controls whether padding is added at the 5' end
(``5``), the 3' end (``3``), or split across both ends (``None``). Padding can
optionally avoid introducing start or stop codons.

Saving results
--------------

.. code-block:: python

    from kea import build_library, save_library

    results = build_library(sequences, "yeast")
    save_library(results, "optimized_sequences.csv")
