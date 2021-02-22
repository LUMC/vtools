==========
Changelog
==========

.. Newest changes should be on top.

.. This document is user facing. Please word the changes in such a way
.. that users understand how the changes affect the new version.

version 2.0.0-dev
-----------------
Numerous backwards incompatible changes have been made to the vtools-gcoverage
script.

+ vtools-coverage's ``--per-transcript`` flag now exhibits different behavior.
  It now gets the coverage for the entire transcript as supplied by the refflat
  file. The old behavior was to only take the exons in the coding region and
  average them. This behavior is still available with the
  ``--per-transcript-cds-exons`` flag.
+ A ``-s`` / ``--short-column-names`` flag was added to the CLI of
  vtools-gcoverage. It sets shorter column names which may enhance readability
  in ``less`` or other terminal readers.
+ The column order of the vtools-coverage output is no longer in alphabetic
  order but in a more human-readable order: means, medians,
  at least percentages depth, at least percentages genome quality.
+ The percentage decimals in the output of vtools-gcoverage are rounded
  to the nearest two digits for better readability.
+ vtools-gcoverage no longer requires a ``-I`` flag to denote the input
  GVCF(s).
+ Multiple single sample GVCFs can now be used as input for vtools-gcoverage.
  The output statistics will be the average over the input GVCFs.
+ Fixed a bug in vtools-gcoverage where all bases from all variants in a
  certain region were considered for coverage statistics, while the first
  and last variant could have bases extruding from the region of interest.
  This means that results processed with this newer version of vtools-gcoverage
  will be slightly different for most of the data.