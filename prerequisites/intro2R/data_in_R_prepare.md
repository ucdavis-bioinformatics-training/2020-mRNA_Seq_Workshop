# Create a new RStudio project

Open RStudio and create a new project, for more info see (Using-Projects)[https://support.rstudio.com/hc/en-us/articles/200526207-Using-Projects]

* File > New Project > New Directory > New Project (name the new directory, Ex. Adv_Mapping_Comparison) and check "use packrat with this project", or "use renv with this project" if your using the devel version.

Learn more about (renv)[https://rstudio.github.io/renv/articles/renv.html]

Learn more about (packrat)[https://rstudio.github.io/packrat/]

## Install Needed Packages

Set some options and make sure the packages 'knitr', 'tidyverse', 'reshape2', and 'gridExtra' are installed (if not install it), and then load

In the R console run the following commands:

```r
if (!any(rownames(installed.packages()) == "knitr")){
  install.packages("knitr")
}
library(knitr)

if (!any(rownames(installed.packages()) == "tidyverse")){
  install.packages("tidyverse")
}
library(tidyverse)

if (!any(rownames(installed.packages()) == "reshape2")){
  install.packages("reshape2")
}
library(reshape2)

if (!any(rownames(installed.packages()) == "gridExtra")){
  install.packages("gridExtra")
}
library(gridExtra)
```

Learn more about the [tidyverse](https://www.tidyverse.org).

## Open a new R Notebook

An [R notebook](https://rmarkdown.rstudio.com/r_notebooks.html) is an R Markdown document with chunks that can be executed independently and interactively, with output visible immediately beneath the input. This is a part of literate programming, there the 'code', description of the 'code' and output are all together in one document.

* File -> New File -> R Notebook
* Save the Notebook (Ex. test)

## R Markdown

This is an R Markdown document. Markdown is a simple formatting syntax for authoring HTML, PDF, and MS Word documents. See [this page](http://rmarkdown.rstudio.com) more details on using R Markdown.

When you click the **preview** or **Knit** button a document will be generated that includes both content as well as the output of any embedded R code chunks within the document. You can embed R code and plots in chunks like this:

<pre><code>```{r chunk_name}
print('hello world!')
```</code></pre>

Review the [R Markdown page](http://rmarkdown.rstudio.com) and R Markdown [cheat sheets](../../base/cheatSheetIndex).

Try 'knitting' to html, pdf, and doc as well as previewing the notebook. Open the resulting documents.

Try executing the code chunks in the R Notebook.

## Download the data file for the workshop document and preview/open it

This is the stats file generated after running samtools stats on a bam file produced from running BWA MEM.

In the R console run the following command.
```r
download.file("https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2020-mRNA_Seq_Workshop/master/prerequisites/intro2R/Data_in_R_files/bwa_mem_Stats.log", "bwa_mem_Stats.log")
```

## Download the template Markdown workshop document and open it

In the R console run the following command
```r
download.file("https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2020-mRNA_Seq_Workshop/master/prerequisites/intro2R/data_in_R.Rmd", "data_in_R.Rmd")
```

## Edit the file YAML portion

The top YAML (YAML ain't markup language) portion of the doc tells RStudio how to parse the document.

<pre><code>---
title: "Data_in_R"
author: your_name
date: current_date
output:
    html_notebook: default
    html_document: default
---</code></pre>

## What are we going to do?

We will recreate some of the plots generated with plot-bamstats on the same file

You can view the output of plot-bamstats -> [bwa_mem_stats.html](Data_in_R_files/bwa_mem_Stats/bwa_mem_Stats.html)
