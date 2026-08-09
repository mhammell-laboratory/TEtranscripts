# Reproducible real count-table fixtures

These 120-feature fixtures exercise compatibility across realistic expression
depths and effect sizes without committing full public datasets.

| Fixture | Public source | Selected libraries | Source SHA-256 |
| --- | --- | --- | --- |
| `pasilla_subsample.cntTable` | [DESeq2 `pasilla_gene_counts.tsv.gz`](https://github.com/thelovelab/DESeq2/blob/devel/inst/extdata/pasilla_gene_counts.tsv.gz) | `treated1`–`treated3` versus `untreated1`–`untreated3` | `4b5570008936786da365d110df1c83c8cffc646d6c667d6b15ecbfb3d9ee735d` |
| `airway_subsample.cntTable` | [Himes airway scaled counts](https://bdsr.stephenturner.us/data/airway_scaledcounts.csv), from [GSE52778](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52778) | dexamethasone-treated `SRR1039509`, `SRR1039513`, `SRR1039517` versus paired untreated `SRR1039508`, `SRR1039512`, `SRR1039516` | `4d34df15affa22fa5d5539874e148378ea18db1c9d731d175ef1b9f73cbc2927` |

`prepare_real_fixtures.py` verifies each downloaded checksum, retains rows that
pass TEtranscripts' default `--minread 1` filter, and stratifies them by four
total-count depth quartiles and three pseudocount log2-fold-change classes
(`down < -1.5`, `stable`, `up > 1.5`). It chooses ten rows per stratum using a
SHA-256 ordering of dataset and feature identifiers. If a stratum is smaller,
the script fills to exactly 120 rows from the remaining stable hash order.

Regenerate or verify the committed byte-for-byte output with:

```console
python tests/compatibility/prepare_real_fixtures.py
python tests/compatibility/prepare_real_fixtures.py --check
```

Regeneration requires network access; CI consumes the committed fixtures and
therefore does not depend on either upstream host being available.
