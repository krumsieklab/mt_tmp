MetaboTools conventions

This document contains conventions on the internal data structures of MetaboTools. It's used as a central reference document.

# SummarizedExperiment

* We use SummarizedExperiment the way it is intended by the authors, i.e. with samples=columns, and rows=variables. Be aware of this if you take the data out of it for non-mt code.
* rowData must have a $name field, which is a printable version of the metabolite/protein’s name 

[incomplete]