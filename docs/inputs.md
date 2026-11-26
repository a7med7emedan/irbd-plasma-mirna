# Input formats

All tables are tab-separated text files with a header row. The folder `examples/data` holds a
synthetic dataset in these formats.

## Count matrix

One row per mature miRNA and one column per library.

| Column | Content |
|---|---|
| first column | miRNA name, for example `hsa-miR-451a`. Names must be unique. |
| other columns | counts per library. The header is the library identifier. |

Rules:

1. Counts must be complete and not negative.
2. Counts may be non-integer, as produced by miRMaster for reads shared between miRNAs. The
   pipeline rounds them.
3. Use miRBase release 22 names, so that the brain-enriched panel and the 14q32 table match.
4. The haemolysis index needs `hsa-miR-451a` and `hsa-miR-23a-3p`.

## Sample sheet

One row per sequenced library.

| Column | Required | Content |
|---|---|---|
| `library_id` | yes | Library identifier, as in the count matrix header |
| `group` | yes | One of the four group labels of the configuration |
| `age` | yes | Age in years |
| `sex` | yes | Two levels, for example `female` and `male` |
| `batch` | yes | Sequencing batch label |
| `qc_pass` | yes | `true` when the library passed quality control |
| `converter` | no | `true` for iRBD individuals who later converted |

Default group labels are `iRBD`, `de_novo_PD`, `HC_high_PRS` and `HC_low_PRS`. Change them under
`groups` in the configuration when your labels differ.

## Clinical table (optional)

One row per analysed library, with `library_id` and one column per clinical variable. Name the
variables under `clinical.continuous` and `clinical.categorical` in the configuration. Categorical
variables need two levels. Empty cells are missing values.

## Resource tables

Two tables ship with the package and are used by default.

1. `brain_enriched_mirnas.tsv`: the 30 predefined brain-enriched miRNAs.
2. `mirbase_v22_chr14_DLK1_DIO3.tsv`: miRNA genes of the DLK1-DIO3 domain with their mature miRNAs,
   derived from miRBase release 22 (GRCh38).

To rebuild the second table from a miRBase GFF3 file:

```bash
irbd-mirna build-cluster-resource --gff hsa.gff3 --out my_domain.tsv
```

Then set `resources.cluster_14q32` to the new file.
