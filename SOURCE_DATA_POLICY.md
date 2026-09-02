# Source-data policy

This release contains no raw or processed input matrices and no identifiable clinical data. It includes only de-identified, panel-level derived tables needed to audit the numerical values displayed in the manuscript figures and supplementary figures.

The public input resources, accession identifiers, modalities and source URLs are listed in `metadata/data_accessions.tsv`. Users must obtain raw or processed inputs from their originating repositories and comply with the respective access conditions.

The included `source_data/` tables are frozen outputs from the documented computational workflow. They are sufficient for figure-level numerical audit, but they do not reconstruct every intermediate object required by the historical provenance scripts. In particular, this archive does not redistribute AnnData objects, source count matrices, protected participant-level clinical information or vendor data.

OEP001105 is hosted through BioSino/NODE; it is not a ProteomeXchange accession. The release manifest identifies the versioned materials included in this archive.
