# Reference required for exact reproduction

Place the custom small-RNA reference here as:

```text
human_hg38_smrna.fa
```

Its required SHA256 is:

```text
b0571c3a903633695f4a27a48d59ecf0b528852389f6e4fee91571e79db731d3
```

This is the reference used to produce the existing HeLa folder-14 results.
It contains tRNA plus the small-RNA/oligo sequences used by the original
alignment, so it cannot be replaced by a generic hg38 tRNA FASTA. The pipeline
checks the hash and builds the Bowtie 1 and samtools indexes automatically.
