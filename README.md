## NOTE
Extensive alpha testing has been conducted, though this software is in a beta state. Errors are expected, often rerunning without changing parameters is sufficient to resume appropriately. 
Kindly raise git issues for errors - if you can find the bug, even better!
Documentation is currently in the works.

<br />

## PURPOSE
The most common gene cluster detection algorithms focus on canonical “core”
biosynthetic functions many gene clusters encode, while overlooking uncommon or
unknown cluster classes. These overlooked clusters are a potential source of
novel natural products and comprise an untold portion of overall gene cluster
repertoires. Unbiased, function-agnostic detection algorithms therefore provide
an opportunity to reveal novel classes of gene clusters and more broadly define
genome organization. *CLOCI* (Co-occurrence Locus and Orthologous Cluster
Identifier) is an algorithm that identifies gene clusters using multiple
proxies of selection for coordinated gene evolution. In the process, *CLOCI*
circumscribes loci into homologous locus groups, which is an extension of
orthogroups to the locus-level. Our approach generalizes gene cluster detection and gene cluster family circumscription, improves detection of multiple known functional classes, and unveils noncanonical gene clusters. *CLOCI* is suitable for genome-enabled specialized metabolite mining, and presents an easily tunable approach for delineating gene cluster families and homologous loci.

<br />

## INSTALL
Please create a conda environment and manually install some dependencies
```bash
conda create -n cloci mycotools graph-tool python pip
```

Then install cloci into the environment
```bash
conda activate cloci
python3 -m pip install cloci
```

A conda package will be available in the future.

<br />

## USE

### Input dataset
*CLOCI* inputs a tab-delimitted file of genome metadata with the following columns: 

```
#genus	species	strain	assembly_path	gffpath
```

or a preassembled [MycotoolsDB](https://github.com/xonq/mycotools). It is important
to adequately sample a cluster's distribution to detect it. I thus generally 
recommend implementing *CLOCI* at least at the subphylum-level. This varies
depending on the lineage's rate of microsynteny decay and the phylogenetic distance 
with which horizontal transfer occurs. 

*CLOCI* fundamentally relies on reconstructing an microsynteny phylogeny that
accurately depicts divergence in gene order between genomes. While *CLOCI*
attempts to automatically detect near single-copy gene families for
reconstructing this tree, it is recommended to explicitly input these
near single-copy genes using the `-f` argument referencing a file of reference genes,
separated by lines. Ideally, these same genes would be used to reconstruct a
phylogenomic tree and the microsynteny topology will be constrained to this
reference phylogenomic tree via the `-c` argument in conjunction with `-r`
for selecting genomes to derive the outgroup branch from.

### Hyperparameters
*CLOCI* default parameters have been tuned for our initial dataset on ~2,250
fungi across the kingdom. These should suffice for circumscribing homologous
loci in most analyses, though the gene cluster family filtering parameters
are ideally determined referencing known clusters from your particular dataset. 
By default, thresholds for all proxies of coordinated gene evolution are set to 0. 
These thresholds will vary for the type of clusters of interest and the
lineage. I recommend compiling a dataset of known cluster reference genes,
running *CLOCI*, identifying those genes in the output, determining the
values for the reference cluster proxies, and then implementing the thresholds.

There are numerous hyperparameters that will drastically affect output quality. 
I suspect our pilot study reached a local maximum in terms of output quality, 
though a global maximum perhaps lies with further hyperparameter tuning. 

### Example
Extract a MycotoolsDB of Agaricomycotina
```bash
mtdb e -l Agaricomycotina > agaricomycotina.mtdb
```

Run *CLOCI* rooting upon the MRCA of two inputted genomes
```bash
cloci -d agaricomycotina.mtdb --root "<OME1>,<OME2>"
```

Resume a *CLOCI* run, i.e. to add proxy thresholds or resume following error
```bash
cloci -d agaricomycotina.mtdb -r <ROOT_OME> -o <PREVIOUS_DIR>
```

<br /><br />


## ON THE ALGORITHM
### Pipeline
![*CLOCI*](https://gitlab.com/xonq/cloci/-/raw/master/etc/pipeline.png)

### Recovery of 68 reference clusters
![Recovery of 68 reference clusters](https://gitlab.com/xonq/cloci/-/raw/master/etc/recovery.png)

### Boundary assessment of 33 reference clusters
![Boundary assessment of 33 reference
clusters](https://gitlab.com/xonq/cloci/-/raw/master/etc/boundaries.png)

### Common Errors
#### Memory error
```
OSError: [Errno 12] Cannot allocate memory
```

Simply resume specifying the run output directory in your command via 
`-o <PREVIOUS_OUTPUT_DIR>`


#### Single-copy gene detection
```
ERROR: could not detect 10 genes present in all genomes with median 2 copy number and less than 2 copy number standard deviation. Manually input focal homology groups.
```

Near single-copy genes were not automatically determined from the dataset. It is recommended to manually input a list of focal homology groups/genes via `-f`.
