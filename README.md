## NOTE
Extensive alpha testing has been conducted, though this software is in a beta state. Errors are expected, often rerunning without changing parameters is sufficient to resume appropriately. 
Kindly raise git issues for errors - if you can find the bug, even better! This software will have longterm maintenance, but maintenance is currently not the priority. Documentation is currently lacking, please be patient as I build out the repository and wiki.

Please also note that this software interfaces with the comparative genomics software suite, Mycotools. I am hopeful you will find Mycotools useful. Please find its manuscript and the [associated repository](https://gitlab.com/xonq/mycotools).

<br />

## PURPOSE
The most common gene cluster detection algorithms focus on canonical “core” biosynthetic functions many gene clusters encode, while overlooking uncommon or unknown cluster classes. These overlooked clusters are a potential source of novel natural products and comprise an untold portion of overall gene cluster repertoires. Unbiased, function-agnostic detection algorithms therefore provide an opportunity to reveal novel classes of gene clusters and more broadly define genome organization. *CLOCI* (Co-occurrence Locus and Orthologous Cluster Identifier) is an algorithm that identifies gene clusters using multiple proxies of selection for coordinated gene evolution. Our approach generalizes gene cluster detection and gene cluster family circumscription, improves detection of multiple known functional classes, and unveils noncanonical gene clusters. CLOCI is suitable for genome-enabled specialized metabolite mining, and presents an easily tunable approach for delineating gene cluster families and homologous loci.

![Recovery of 68 reference clusters](https://gitlab.com/xonq/cloci/-/raw/master/etc/recovery.png)

![Boundary assessment of 33 reference
clusters](https://gitlab.com/xonq/cloci/-/raw/master/etc/boundaries.png)


<br /><br />

## APPROACH
![CLOCI](https://gitlab.com/xonq/cloci/-/raw/master/etc/pipeline.png)

