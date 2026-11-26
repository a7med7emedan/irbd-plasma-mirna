# Workflow

![Workflow of the pipeline, from the inputs to the key numbers](img/workflow.png)

The figure shows every step of `irbd-mirna run`, in the order of `pipeline.run`. Labels and numbers
come from [methods.md](methods.md), `pipeline.py`, `de.py` and `r/deseq2_functions.R`. The defaults
shown are those of `config.py`. The two small plots in the figure are schematic and contain no data.

The figure is drawn by `tools/make_doc_figures.py`, together with `img/how-it-works.png`. Run
`make doc-figures` after a change to the pipeline.

<details>
<summary>Text version of the diagram</summary>

```mermaid
flowchart TB
    subgraph S1 ["1. Inputs"]
        direction LR
        I1["Count matrix<br/>miRNA by library<br/>miRMaster 2.0 counts"]
        I2["Sample sheet<br/>library_id, group, age, sex,<br/>batch, qc_pass, converter"]
        I3["Clinical table<br/>optional"]
        I4["brain_enriched_mirnas.tsv<br/>30 predefined miRNAs"]
        I5["mirbase_v22_chr14_DLK1_DIO3.tsv<br/>miRBase v22, GRCh38"]
        I1 ~~~ I2 ~~~ I3 ~~~ I4 ~~~ I5
    end

    subgraph S2 ["2. Input checks and quality control, Python, io.py"]
        direction LR
        Q1["Input checks<br/>required columns, known group labels,<br/>numeric and unique count columns"]
        Q2["Quality-control filter<br/>libraries with qc_pass"]
        Q3["Rounding<br/>non-integer counts to integers"]
        Q4["Detection summary and read depth<br/>reads per million mapped miRNA reads"]
        Q1 --> Q2 --> Q3 --> Q4
    end

    subgraph S3 ["3. Differential expression, R, DESeq2, run_deseq2.R"]
        direction TB
        subgraph S3A ["Four comparisons, each fitted on its own two groups"]
            direction LR
            C1["High-PRS vs low-PRS controls<br/>pre-filter: edgeR filterByExpr"]
            C2["iRBD vs controls<br/>pre-filter: edgeR filterByExpr"]
            C3["iRBD vs de novo PD<br/>pre-filter: edgeR filterByExpr"]
            C4["de novo PD vs controls<br/>pre-filter: at least 10 counts<br/>in at least 3 libraries"]
            C1 ~~~ C2 ~~~ C3 ~~~ C4
        end
        subgraph S3B ["Fitted for each comparison"]
            direction LR
            D1["Model ~ age + sex + group<br/>age centred and scaled<br/>no batch term"]
            D2["Median-of-ratios size factors<br/>dispersions shrunk<br/>towards a fitted trend"]
            D3["Wald test"]
            D4["Independent filtering, alpha 0.1<br/>Benjamini-Hochberg adjustment"]
            D5["apeglm shrinkage<br/>of log2 fold changes"]
            D1 --> D2 --> D3 --> D4 --> D5
        end
        D6["Also written: DESeq2 size factors over all libraries and within iRBD"]
        S3A --> S3B
        S3B ~~~ D6
    end

    subgraph S4 ["4. Dysregulation rule, Python, dysregulation.py"]
        R1["Dysregulated when FDR below 0.05<br/>and absolute shrunken log2 fold change at least 0.5"]
    end

    subgraph S5 ["5. Downstream analyses, Python"]
        direction TB
        subgraph S5A [" "]
            direction LR
            A1["14q32 cluster<br/>by genomic position, forward strand of<br/>chromosome 14 between DLK1 and DIO3<br/>tested and dysregulated per comparison<br/>hsa-mir-1247 reported separately"]
            A2["Brain-enriched panel<br/>status in the three main comparisons<br/>shared candidates of both iRBD comparisons<br/>one-sided Fisher's exact test"]
            A3["Haemolysis index and cluster score<br/>log2 of miR-451a over miR-23a-3p<br/>two-sided Mann-Whitney U, Spearman<br/>regressions: group, plus index, plus age and sex"]
            A1 ~~~ A2 ~~~ A3
        end
        subgraph S5B [" "]
            direction LR
            A4["Converters within iRBD<br/>detected in at least 20% of iRBD libraries<br/>two-sided Mann-Whitney U<br/>200 permuted datasets"]
            A5["Batch analysis<br/>chi-square test and Cramer's V<br/>principal components, Kruskal-Wallis<br/>per-miRNA F test, Benjamini-Hochberg"]
            A6["Clinical comparisons<br/>Welch's t test, Fisher's exact test<br/>Benjamini-Hochberg, missing data"]
            A4 ~~~ A5 ~~~ A6
        end
        S5A ~~~ S5B
    end

    subgraph S6 ["6. Figures, Python, figures.py"]
        direction LR
        F1["Figure 1<br/>volcano plots of the<br/>three main comparisons"]
        F2["Figure 2<br/>brain-enriched panel"]
        F3["Supplementary Figure S1<br/>high-PRS vs low-PRS controls"]
        F4["Supplementary Figure S2<br/>converters vs non-converters"]
        F5["Supplementary Figure S3<br/>batch principal components"]
        F1 ~~~ F2 ~~~ F3 ~~~ F4 ~~~ F5
    end

    subgraph S7 ["7. Key numbers, Python, report.py"]
        direction LR
        K1["summary/key_numbers.tsv<br/>summary/run_info.yaml"]
        K2(["Compare with<br/>config/manuscript_values.yaml"])
        K3["match, differs or not reported<br/>--strict stops on differs"]
        K1 --> K2 --> K3
    end

    subgraph S8 ["Reference image"]
        direction LR
        V1["R 4.6, Bioconductor 3.23<br/>DESeq2 1.52, apeglm 1.34, edgeR"]
        V2["Python 3.10.12<br/>SciPy 1.15.3, statsmodels 0.14.6"]
        V1 ~~~ V2
    end

    S1 --> S2
    S2 --> S3
    S3 --> S4
    S4 --> S5
    S5 --> S6
    S6 --> S7
    S7 ~~~ S8

    classDef rstep fill:#dbeafe,stroke:#1e40af,color:#0f172a
    classDef pystep fill:#fef3c7,stroke:#92400e,color:#0f172a
    classDef input fill:#f1f5f9,stroke:#475569,color:#0f172a
    classDef version fill:#ffffff,stroke:#475569,color:#0f172a,stroke-dasharray:4 3
    class I1,I2,I3,I4,I5 input
    class C1,C2,C3,C4,D1,D2,D3,D4,D5,D6 rstep
    class Q1,Q2,Q3,Q4,R1,A1,A2,A3,A4,A5,A6,F1,F2,F3,F4,F5,K1,K2,K3 pystep
    class V1,V2 version
    style S5A fill:none,stroke:none
    style S5B fill:none,stroke:none
```

</details>

## Legend

1. In the figure, blue tags mark steps in R and amber tags steps in Python. Red frames mark the
   dysregulation rule and the key numbers. Dashed cards are optional. In the text version, blue
   boxes run in R, amber boxes run in Python and grey boxes are inputs.
2. Arrows between stages follow the order of `pipeline.run`. The 14q32 cluster is counted before
   the brain-enriched panel, because the panel test also uses a background without the cluster.
3. Stage 5 also reads the checked inputs of stage 2. The haemolysis index and the batch analysis use
   reads per million mapped miRNA reads. The converter comparison uses the DESeq2 size factors
   within iRBD. The clinical comparisons use the clinical table.
4. Controls are the high-PRS and low-PRS controls pooled. Fold changes are the first-named group
   relative to the second-named group.
5. The batch analysis keeps high-PRS and low-PRS controls apart and uses the leading three
   principal components.
6. The permutation check uses the seed 20260830. The comparison with the manuscript values is
   switched off for the synthetic example data.
