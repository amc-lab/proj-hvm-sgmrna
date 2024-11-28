## Summary

This repository contains the Nextflow pipeline and run script to process the amplicons sequencing data for: Mears et. al (2023) __Emergence of new subgenomic mRNAs in SARS-CoV-2__

### Dependencies

This pipeline was run using:
- Nextflow v23.10.0
- Singularity v3.6.4

## Abstract

Coronaviruses express their structural and accessory genes via a set of subgenomic RNAs, whose synthesis is directed by Transcription Regulatory Sequences (TRS) in the 5′ genomic leader and upstream of each body open reading frame. In SARS-CoV-2, the TRS has the consensus AAACGAAC; upon searching for emergence of this motif in the global SARS-CoV-2 sequences, we find that it evolves frequently, especially in the 3′ end of the genome. We show well-supported examples upstream of the Spike gene – within the nsp16 coding region of ORF1b – which is expressed during human infection, and upstream of the canonical Envelope gene TRS, both of which have evolved convergently in multiple lineages. The most frequent neo-TRS is within the coding region of the Nucleocapsid gene, and is present in virtually all viruses from the B.1.1 lineage, including the variants of concern Alpha, Gamma, Omicron and descendants thereof. Here, we demonstrate that this TRS leads to the expression of a novel subgenomic mRNA encoding a truncated C-terminal portion of Nucleocapsid, which is an antagonist of type I interferon production and contributes to viral fitness during infection. We observe distinct phenotypes when the Nucleocapsid coding sequence is mutated compared to when the TRS alone is ablated. Our findings demonstrate that SARS-CoV-2 is undergoing evolutionary changes at the functional RNA level in addition to the amino acid level.