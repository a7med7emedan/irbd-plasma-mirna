# Troubleshooting

| Message | Cause | What to do |
|---|---|---|
| `'Rscript' was not found` | R is not installed or not on the path | Install R with DESeq2, edgeR and apeglm, set `rscript` in the configuration, or use the Docker image |
| `The DESeq2 step failed` | R stopped with an error | Read the lines starting with `R |` in the log. Common causes are a covariate confounded with group and a group without libraries |
| `The sample sheet lacks these columns` | A required column is missing | Add the column, see [inputs.md](inputs.md) |
| `The sample sheet has unknown group labels` | Labels differ from the configuration | Change the labels, or set them under `groups` |
| `libraries of the sample sheet are not in the count table` | Identifiers differ between the two tables | Make `library_id` equal to the count matrix header |
| `Unknown configuration key` | A key is misspelled or not supported | Compare with [usage.md](usage.md) |
| `hsa-miR-451a is not in the count table` | The haemolysis index needs this miRNA | Use miRBase release 22 names in the count matrix |
| `Converter comparison skipped` | Fewer than two libraries in one group | Check the `converter` column, or ignore when the step is not needed |
| `Batch analysis skipped` | Only one batch | No action needed |
| `key numbers differ from the manuscript values` | `--strict` found differences | Open `summary/key_numbers.tsv` and check the inputs and the software versions in `summary/run_info.yaml` |
| `Python 3.10 or newer is required` | The Python version is too old | Use the conda environment or the Docker image |
| Note about `fitType='parametric'` in the log | DESeq2 switched to a local dispersion fit | No action needed. This is expected on small or synthetic datasets |

## Getting help

Open an issue in the repository. Include the command, the last 30 log lines and the file
`summary/run_info.yaml`. Do not attach participant data.
