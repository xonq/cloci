## Repository setup
- [ ] Begin USAGE wiki
- [ ] Begin INSTALL wiki
- [ ] Begin Output wiki
- [x] Build pypi package
- [ ] Build conda package
- [x] Update intro page
- [ ] Add issue format
- [x] Mirror with GitHub

## Algorithm improvements
- [ ] Merge loci extensions that overlap after the second round
- [ ] Option to merge adjacent reported loci that belong to the same HLG on
  output based on a gene intermediance
- [x] Fix rooting system for nodes that overlap current root/ambiguous
- [ ] Multiprocess null generation
- [ ] Lineage-based thresholds
- [ ] Standardize code annotation
- [ ] Account for tandem dup HGs
- [ ] Tune hyperparameters for global maximum quality output
	- [ ] High percentage of 1 domain/reference cluster
- [ ] Convert sliding window to nucleotide-based slide
- [ ] Account for alternate splicing appropriately (gene-wise v RNA-wise)
- [ ] Implement random sampling probability HGps and HGxs
- [ ] Account for duplicates of HGs on edge of cluster
- [ ] Add Pfam extraction feature
- [ ] Add GO extraction feature
- [ ] Output networks as interactable files
- [ ] Phylogeny-based method of GCL calculation
- [ ] Add a force skip to filtering that reads in old runs and converts
- [ ] Output version in the log, have a version output statement
- [ ] aPDS that removes the signal of HLG loss
- [ ] rPDS calculation with GCL threshold in cloci2stats

## Quality of life
- [ ] pruning discrepancies from microsynteny tree and input database and
  removing from analysis
- [ ] common visualization outputs
