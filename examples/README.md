# Example data

The files in `data/` are synthetic and represent no participant. They show the input formats and
let the pipeline run end to end.

| File | Content |
|---|---|
| `data/counts.tsv` | 413 miRNAs by 52 libraries |
| `data/samples.tsv` | Sample sheet with four groups, 48 libraries passing quality control |
| `data/clinical.tsv` | Three clinical variables with some missing values |

The dataset has planted effects, so that every step has something to find:

1. The forward-strand miRNAs of the 14q32 cluster are lower in iRBD.
2. Thirty background miRNAs and three panel miRNAs are higher in iRBD.
3. Three background miRNAs are higher in de novo PD.
4. High-PRS and low-PRS controls do not differ.

Rebuild the files with:

```bash
make example-data
```

Run the pipeline on them with:

```bash
make run
```
