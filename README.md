## NOTE
While extensive alpha testing has been conducted, this software is in a beta state, and errors are expected. Kindly report these issues and remain patient for fixes - if you can find the bug, even better! This software will have longterm maintenance, but maintenance is currently not the priority. Documentation is currently lacking, please be patient as I build out the repository and wiki.

Please also note that this software interfaces with the comparative genomics software suite, Mycotools, which entails its own setup and learning curve. I am hopeful you will find Mycotools useful. Please find its manuscript and the [associated repository](https://gitlab.com/xonq/mycotools).

Cheers,
Zachary Konkel

## PURPOSE
Identifying loci comprised of cooperative biological functions, or gene clusters, can pinpoint biosynthetic pathways for drug discovery and provide insight into genomes’ ecological functions. The first widely-adapted gene cluster prediction algorithms are largely function-centric, which screen for modeled, common biosynthetic gene cluster functions in new genomes. While these approaches have enabled a genome-guided era of drug discovery, function-centric detection cannot de novo infer noncanonical gene cluster classes which lack modeled functions. Noncanonical clusters, such as the psilocybin gene cluster, thus present a significant oversight in drug discovery research. Additionally, the function-centric paradigm is inherently biased toward well-studied organisms’ biosynthetic gene cluster classes and overlook other ecologically-relevant types of gene clusters. We present Homology Groups to Clusters (*Cloci*) as an evolution-based, generalized gene cluster detection pipeline that can infer noncanonical gene cluster classes de novo. *Cloci* identifies families of related gene clusters by aggregating loci that are syntenic across an unexpected sample of genomes and filtering groups of homologous loci using measurements of selection for gene clustering. We benchmarked *Cloci* against contemporary cluster detection algorithms and show that *Cloci* has the highest recall and cluster boundary inference, alongside functionally-significant cluster family aggregation. Our benchmarked dataset includes catabolic, nutrient assimilation, and biosynthetic gene clusters, which demonstrates that *Cloci* is poised for adoption in comparative ecological analyses and generalized gene cluster detection. *Cloci* presents a paradigm shift in gene cluster detection from function-centric to incorporating function-agnostic detection as a viable orthogonal approach to comparative genomics and drug discovery research.

<br /><br />

## APPROACH
![Cloci](https://gitlab.com/xonq/orthocluster/-/raw/master/etc/pipeline.png)

