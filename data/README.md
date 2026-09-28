# Input data

Nothing in this folder is committed (see `.gitignore`). Download the files below
and place them at these paths; the pipeline scripts use them by default.

| Path | What it is | Where to get it |
|---|---|---|
| `kinome_atlas/ST2.xlsx` | Position-specific matrices of 303 Ser/Thr kinases (sheet `ser_thr_all_norm_scaled_matrice`) | Johnson *et al.*, *Nature* (2023), Supplementary Table 2 — https://doi.org/10.1038/s41586-022-05575-3 |
| `kinome_atlas/ST3.csv` | Reference human phosphoproteome with published kinase percentiles (save Supplementary Table 3 as CSV) | same paper, Supplementary Table 3 |
| `mouse_atlas/ST2.xlsx` | Mouse phosphosites across 41 tissues (demo input; sheet `Supplementary Table 2 - Organs`) | Giansanti *et al.*, *Nat. Methods* (2022), Supplementary Table 2 — https://doi.org/10.1038/s41592-022-01526-y (raw data: PRIDE PXD030983) |
| `uniprot/mouse.tsv` | Mouse proteome | UniProt → *Mus musculus* (taxon 10090), download TSV with columns Entry, Reviewed, Gene Names (primary), Sequence |
| `uniprot/human.tsv` | Human proteome | UniProt → *Homo sapiens* (taxon 9606), same columns |

Example UniProt queries (REST):

```
https://rest.uniprot.org/uniprotkb/stream?format=tsv&fields=accession,reviewed,gene_primary,sequence&query=organism_id:10090
https://rest.uniprot.org/uniprotkb/stream?format=tsv&fields=accession,reviewed,gene_primary,sequence&query=organism_id:9606
```

For your own experiment, replace the mouse atlas with MaxQuant's
`Phospho (STY)Sites.txt` (step 1) and a limma table per comparison (step 4).
