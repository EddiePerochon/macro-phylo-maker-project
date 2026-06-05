# Project To-Do List: MacroPhyloMaker

This document tracks planned features, improvements, and maintenance tasks.

---

## High Priority

- [ ] Write unit tests for all exported functions (≥80% coverage)
- [ ] Make sure AntWiki taxonomy to TACT template works
- [ ] Correctly handle paraphyletic genera in TACT templates
  - [ ] Oxyopomyrmex: 
        Nested in Goniomma, 2 species in backbone but 14 total; 
        - [ ] Can graft separately using TACT, 
        - [ ] prune from backbone, 
        - [X] and replace back onto final tree sister to Goniomma blanci
  - [ ] Aphaenogaster: 
        Tropical Aphaenogaster clade can be called Deromyrma in both tree and csv;
        - [ ] This requires figuring out which species missing from the phylogeny belong to Deromyrma
        hen Aphaenogaster s. str. is paraphyletic with respect to Messor;
        - [ ] But in Juve et al. 2025 it is not. So this can be fixed by extracting both Aphaenogaster and Messor from that tree
  - [ ] Rossomyrmex:
        Nested in Cataglyphis.
        - [ ] Can graft separately using TACT, 
        - [ ] prune from backbone, 
        - [X] and replace back onto final tree?
- [ ] Rectify duplicate Strumigenys emmae
- [ ] Rectify duplicate Mystrium mysticum
- [ ] Rectify duplicate Mystrium rogeri

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

---