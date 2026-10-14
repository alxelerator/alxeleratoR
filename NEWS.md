# alxeleratoR 0.8.4.9200 (2026-10-14)

* Added automatic italic formatting of species and genus names. Override can turn this off. 
  It respects the other rank identifiers at the end of i.e. GenusSpecies names [Family] etc...

# alxeleratoR 0.8.4.9000 (dev)

* Fixed rarefy_custom_depth(): rewrote loop using lapply and matrix cbind/rbind to
  avoid S4 scoping issue where iterative merge_phyloseq inside a function only returned
  the last sample. Function is now exported and working.
* Added ps_plot_sankey(): Sankey diagram of taxonomic composition using the sankey
  package (Csardi & Weiner). Features: 100% column scaling so all rank axes have
  equal visual height; top-N filter defaults to the finest rank (to_rank) to keep
  intermediate bars manageable; abundance-sorted stacking (largest at bottom);
  Other bundling for non-top taxa; SVG and PNG output via file= argument.
  Replaces a planned ggalluvial-based implementation.

# alxeleratoR 0.8.3.000 (2026-04-22)

* Added functionality to ps_plot_bar to add "other" taxon category in top-x plots

# alxeleratoR 0.8.2.000 (2025-09-03)

* Fixed a bug when using log scales in sample_sum_plot(). Thanks to Tirza for reporting and fixing this bug.
* Fixed many roxygen2 function documentation issues.
* Additional exports of functions added.

# alxeleratoR 0.8.1.000 (2025-04-04)

* First complete build of the package and streamlining documentation.

# alxeleratoR 0.8.0.000 (2025-01-11)

* Initial R package build from existing alxR helper-functions bundle.

# alxeleratoR 0.0.1.000 (2025-01-03)

* Initial checkin of the project.
