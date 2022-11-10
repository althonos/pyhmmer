Easel
=====

.. currentmodule:: pyhmmer.easel


.. automodule:: pyhmmer.easel


Data Structures
---------------

Bitfield
^^^^^^^^

.. autoclass:: pyhmmer.easel.Bitfield
   :special-members: __init__
   :members:

KeyHash
^^^^^^^

.. autoclass:: pyhmmer.easel.KeyHash
   :special-members: __init__
   :members:


Sequences
---------

Sequence
^^^^^^^^

.. autoclass:: pyhmmer.easel.Sequence
   :members:

TextSequence
^^^^^^^^^^^^

.. autoclass:: pyhmmer.easel.TextSequence(Sequence)
   :special-members: __init__
   :members:

DigitalSequence
^^^^^^^^^^^^^^^

.. autoclass:: pyhmmer.easel.DigitalSequence(Sequence)
   :special-members: __init__
   :members:


Sequence Blocks
---------------

SequenceBlock
^^^^^^^^^^^^^

.. autoclass:: pyhmmer.easel.SequenceBlock
   :members:

TextSequenceBlock
^^^^^^^^^^^^^^^^^

.. autoclass:: pyhmmer.easel.TextSequence(SequenceBlock)
   :special-members: __init__
   :members:

DigitalSequenceBlock
^^^^^^^^^^^^^^^^^^^^

.. autoclass:: pyhmmer.easel.DigitalSequenceBlock(SequenceBlock)
   :special-members: __init__
   :members:


SequenceFile
------------

SequenceFile
^^^^^^^^^^^^

.. autoclass:: pyhmmer.easel.SequenceFile
  :special-members: __init__
  :members:


Alignments
----------

MSA
^^^

.. autoclass:: pyhmmer.easel.MSA
   :members:

TextMSA
^^^^^^^

.. autoclass:: pyhmmer.easel.TextMSA(MSA)
   :special-members: __init__
   :members:

DigitalMSA
^^^^^^^^^^

.. autoclass:: pyhmmer.easel.DigitalMSA(MSA)
   :special-members: __init__
   :members:

MSAFile
^^^^^^^

.. autoclass:: pyhmmer.easel.MSAFile
   :special-members: __init__
   :members:


Linear Algebra
--------------

Vector
^^^^^^

.. autoclass:: pyhmmer.easel.Vector
   :members:

VectorF
^^^^^^^

.. autoclass:: pyhmmer.easel.VectorF
   :special-members: __init__
   :members:

VectorU8
^^^^^^^^

.. autoclass:: pyhmmer.easel.VectorU8
   :special-members: __init__
   :members:

Matrix
^^^^^^

.. autoclass:: pyhmmer.easel.Matrix
   :members:

MatrixF
^^^^^^^

.. autoclass:: pyhmmer.easel.MatrixF
   :special-members: __init__
   :members:

MatrixU8
^^^^^^^^

.. autoclass:: pyhmmer.easel.MatrixU8
   :special-members: __init__
   :members:


Miscellaneous
-------------

Alphabet
^^^^^^^^

.. autoclass:: pyhmmer.easel.Alphabet
   :members:

Randomness
^^^^^^^^^^

.. autoclass:: pyhmmer.easel.Randomness
   :special-members: __init__
   :members:

Sequence / Subsequence Index
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: pyhmmer.easel.SSIReader
   :special-members: __init__
   :members:

.. autoclass:: pyhmmer.easel.SSIWriter
   :special-members: __init__
   :members:
