<b>Analysis code for RNA sequencing analysis of P. aeruginosa growth on mucin</b>
This work is in support of Arif and Hoffman et al. (2025)

This repository contains R Markdown notebooks used to process data, perform statistical analyses, and generate figures for the manuscript. The workflow proceeds step-by-step in numbered .Rmd files, and all code can be run interactively or knitted sequentially.

<b>Repository structure</b>
<code>
code/     # R Markdown files for each analysis step (numbered)
data/     # Input for R Markdown files
docs/     # Metadata and other files
</code>

<b>Requirements</b>
<code>
library(DESeq2)
library(phyloseq)
library(tidyverse)
library(pheatmap)
library(clusterProfiler)
library(here)
</code>

<b>How to Run</b>
To reproduce the full analysis:

Open RStudio (recommended).

Navigate to the code/ directory.

Run or knit the scripts in order:
<code>
01_RSubread.Rmd                # Read alignment and quantification
02_DA-analysis.Rmd             # Differential expression analysis
03_functional-analysis.Rmd     # Over-representation analysis
04_heatmap.Rmd                 # Heatmap visualization of key genes
</code>
Each .Rmd script contains its own narrative and code. Parameters are generally hard-coded within the notebooks for transparency.

This repository contains only example or processed data. Raw sequencing data are available at:
[SRA link]

If you use this code, please cite the associated publication:
Author(s). (2025). Title. *Journal*, Volume(Issue), Pages. DOI

License
This project is released under the ___ License.
Thanks to [lab name, collaborators, funding sources].


