# Project To-Do List: MacroPhyloMaker

This document tracks planned features, improvements, and maintenance tasks.

---

## High Priority

- [X] Count genera in each tree
- [X] Find out which genera are missing from species tree with phylogenetic info
- [ ] Fix issue with grafting dropping tips where species name starts with "sp"
- [ ] Figure out how to deal with Gaolingidris: In the all-genera grafted tree under seed 42 it is sister to Malagidris within a clade of Malagasy-endemic genera. This seems unlikely. However, grafting it outside of these genera will create inconsistency between species-grafted and genus-grafted trees. Similar situation may br the case with Peronomyrmex being sister to Pristomyrmex. Is that realistic?
- [ ] Graft missing genera back onto the species tree with phylogenetic info 
- [ ] Redo Amblyoponinae extraction from Nelsen 2018
- [ ] Write unit tests for all exported functions (≥80% coverage)
- [ ] Make sure AntWiki taxonomy to TACT template works
- [ ] Correctly handle paraphyletic genera in TACT templates
  - [ ] Oxyopomyrmex: Nested in Goniomma, 2 species in backbone but 14 total; 
    - [ ] Can graft separately using TACT, 
    - [ ] prune from backbone, 
    - [X] and replace back onto final tree sister to Goniomma blanci
  - [ ] Aphaenogaster: Tropical Aphaenogaster clade can be called Deromyrma in both tree and csv;
    - [ ] This requires figuring out which species missing from the phylogeny belong to Deromyrma
    - [ ] Then Aphaenogaster s. str. is paraphyletic with A. subterranea+pallida closer to Messor; But in Juve et al. 2025 non-Deromyrma Aphaenogaster is monophyletic. But chronogram does not include this lineage. So this can be fixed by pruning A. subterranea and pallida from Branstetter tree OR extracting this lineage from phylogram of Juve et al. 2025.
  - [ ] Rossomyrmex: Nested in Cataglyphis.
    - [ ] Can graft separately using TACT, 
    - [ ] prune from backbone, 
    - [X] and replace back onto final tree?
  - [ ] Stigmatomma: Polyphyletic. The type species is S. denticulatum and it is a clade with S. scrobiceps, S. silvestrii. Would need to remove S. boltoni, S. pallipes, S. oregonense.
    - [ ] Remove Stigmatomma_boltoni, Stigmatomma_pallipes, Stigmatomma_oregonense and
    - [ ] Graft remaining using TACT? This will basically mean we are making it monophyletic depite knowing it is not.
    - [ ] Alternatively, graft Stigmatomma but also prune Adetomyrma, Xymmer, and Mystrium. Then re-graft those three genera within Amblyopone (leaning towards that) 
- [ ] Consider incorporating helper CLI scripts into package (AntWiki to TACT taxonomy, filling genus-only tips with species) 

## Testing & Validation

- [ ] Add unit tests (testthat)
- [ ] Benchmark performance on large datasets
- [ ] Validate output against published phylogenies
- [ ] Create reproducible example workflows


## Documentation

- [ ] Ensure documentation is complete for all functions
- [ ] Build pkgdown website
- [ ] Add “Getting Started” vignette
- [ ] Add advanced use-case vignette
- [ ] Document assumptions and limitations

## Maintenance

- [ ] Continuous Integration setup (GitHub Actions)
- [ ] CRAN submission checks
- [ ] Code style consistency (lintr/styler)
- [ ] Dependency updates

## Done

- [X] Check whether the extract clade function correctly parses names with "cf" in them (e.g., Dorylinae) 
- [X] Ensure extract clade writes dropped tips log
- [X] Rectified duplicate Pheidole nemoralis based on species-group placement in Salata and Fisher 2020 (10.3897/zookeys.905.39592)
- [X] Move helper scripts into package directory
- [X] Rectify duplicate Strumigenys emmae: Removed from MAMI Malagasy tree
- [X] Rectify duplicate Mystrium mysticum: Corrected MRCA placement for MAMI tree
- [X] Rectify duplicate Mystrium rogeri: Corrected MRCA placement for MAMI tree

---