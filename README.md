## NOTE
Extensive alpha testing has been conducted, though this software is in a beta state. Errors are expected, often rerunning without changing parameters is sufficient to resume appropriately. 
Kindly raise git issues for errors - if you can find the bug, even better! This software will have longterm maintenance, but maintenance is currently not the priority. Documentation is currently lacking, please be patient as I build out the repository and wiki.

Please also note that this software interfaces with the comparative genomics software suite, Mycotools. I am hopeful you will find Mycotools useful. Please find its manuscript and the [associated repository](https://gitlab.com/xonq/mycotools).

<br />

## PURPOSE
The most common gene cluster detection algorithms focus on canonical “core” biosynthetic functions many gene clusters encode, while overlooking uncommon or unknown cluster classes. These overlooked clusters are a potential source of novel natural products and comprise an untold portion of overall gene cluster repertoires. Unbiased, function-agnostic detection algorithms therefore provide an opportunity to reveal novel classes of gene clusters and more broadly define genome organization. *CLOCI* (Co-occurrence Locus and Orthologous Cluster Identifier) is an algorithm that identifies gene clusters using multiple proxies of selection for coordinated gene evolution. Our approach generalizes gene cluster detection and gene cluster family circumscription, improves detection of multiple known functional classes, and unveils noncanonical gene clusters. CLOCI is suitable for genome-enabled specialized metabolite mining, and presents an easily tunable approach for delineating gene cluster families and homologous loci.

<br />

## INSTALL
pypi and conda packages will be built for CLOCI in the coming weeks. Until
then, please clone the repository and install the necessary dependencies
manually.

<br />

## USE

### Input dataset
CLOCI inputs a [MycotoolsDB](https://gitlab.com/xonq/mycotools)
that contains genomes ('omes') of interest. It is important
to adequately sample a cluster's distribution to detect it. I thus generally 
recommend implementing CLOCI at least at the subphylum-level. This varies
depending on the lineage's rate of microsynteny decay and the phylogenetic distance 
with which horizontal transfer occurs. 

### Hyperparameters
CLOCI default parameters have been tuned for our initial dataset on ~2,250
fungi across the kingdom. These should suffice for most analyses. By default,
thresholds for all proxies of coordinated gene evolution are set to 0. 
These thresholds will vary for the type of clusters of interest and the
lineage. I recommend compiling a dataset of known cluster reference genes,
running CLOCI, identifying those genes in the output, determining the
values for the reference cluster proxies, and then implementing the thresholds.

There are numerous hyperparameters that will drastically affect output quality. 
I suspect our pilot study reached a local maximum in terms of output quality, 
though a global maximum perhaps lies with further hyperparameter tuning. 

### Example
Extract a MycotoolsDB of Agaricomycotina
```bash
mtdb e -l Agaricomycotina > agaricomycotina.mtdb
```

Run CLOCI
```bash
cloci -d agaricomycotina.mtdb -r <ROOT_OME>
```

Resume a CLOCI run, i.e. to add proxy thresholds or resume following error
```bash
cloci -d agaricomycotina.mtdb -r <ROOT_OME> -o <PREVIOUS_DIR>
```

<br /><br />


## ON THE ALGORITHM
### Pipeline
![CLOCI](https://gitlab.com/xonq/cloci/-/raw/master/etc/pipeline.png)

<br />

### Recovery of 68 reference clusters

![Recovery of 68 reference clusters](https://gitlab.com/xonq/cloci/-/raw/master/etc/recovery.png)


### Boundary assessment of 33 reference clusters
![Boundary assessment of 33 reference
clusters](https://gitlab.com/xonq/cloci/-/raw/master/etc/boundaries.png)
