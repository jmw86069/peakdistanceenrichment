
# Peak-Distance-Enrichment

This protocol describes the process for running distance-based
enrichment methods used in Langer et al 2019, as shown in
a visual workflow in Figure 2, supplement C.

* Provide "reference sites" that have genome coordinates, and
a category assigned to each site.

   * In this case the reference sites are transcript start sites (TSS)
   of genes tested by RNA-seq.
   * The category column named `"hit"` contains
   `-1` for down-regulated genes,
   `1` for up-regulated genes, and
   `0` for genes not changed.
   * Note that "reference sites" can have any category values
   relevant to the experiment, but must contain at least two
   values in order to test enrichment of one category compared
   to the other. In general, "reference sites" should represent
   the universe of sites detectable by the appropriate experimental
   protocol.

* Define "peaks" that will be tested for association to reference sites,
where each peak has genome coordinates. Note for this workflow, strand
is not considered.

   * In this case there are two sets of peaks: ATAC-seq peaks that
   have **higher** chromatin accessibility upon BAF47 knockdown; and
   ATAC-seq peaks that have **lower** chromatin accessibility upon BAF47
   knockdown.
   * In reality any set of peaks can be used, and multiple sets
   can be tested in series.
   
* Start at distance 0 and iterate through the following steps.
Between each iteration, increase the distance by 5000 and repeat.

   * Expand the set of "peaks" by the distance from the center of each peak.
   After this expansion, all "peaks" will be the same width.
   * Merge any overlapping peaks so that there are no overlapping regions.
   This step is sometimes called "merge", "disjoin", or "reduce".
   * Tabulate the number of reference sites that overlap the
   expanded peak regions, for each category value; also tabulate
   the number of reference sites that do not overlap the expanded
   peak regions, for each category value.

      * Using the example data, this step will produce 6 values:

      1. Number of sites with category `-1` that overlap peaks
      2. Number of sites with category `0` that overlap peaks
      3. Number of sites with category `1` that overlap peaks
      4. Number of sites with category `-1` that do not overlap peaks
      5. Number of sites with category `0` that do not overlap peaks
      6. Number of sites with category `1` that do not overlap peaks

* For each distance threshold, an enrichment test is performed
using one or more "test category values", compared to the other
"background category values".

   * For the example above, one may test `1` for enrichment compared
   to `0, -1`.
   * Alternatively, one may test `-1, 1` for enrichment compared to `0`.

* A hypergeometric enrichment test is performed, which requires
the following four values:

   1. Number of sites that overlap peaks, which have "test category values".
   2. Number of sites that overlap peaks, which have "background category values".
   3. Number of sites that do not overlap peaks, which have "test category values".
   4. Number of sites that do not overlap peaks, which have "background category values".

* In R, the function `phyper()` is used to perform this test,
with the following arguments:

   * `q` - number of sites with category of interest that overlap peaks `(1)`
   * `m` - total number of sites with category of interest `(1 + 3)`
   * `n` - total number of sites that do not have the category of interest `(2 + 4)`
   * `k` - number of sites that overlap peaks `(1 + 2)`


## R analysis workflow

There is an R package in a Github repository: https://github.com/jmw86069/peakdistanceenrichment

The R package can be installed by opening an R session and
running this command:

```R
remotes::install_github("jmw86069/peakdistanceenrichment")
```

The package installs four example files, sufficient to
reproduce Figure 2B of Langer et al 2019, eLife.

* `hg19.sizes.txt - chromosome sizes
* `BAF47KD_steadystate_RNAseq_hits.bed` - reference TSS sites with `-1,0,1` hits
* `BAF47KD_steadystate_higher_accessibility_peaks.bed` - test peaks
* `BAF47KD_steadystate_lower_accessibility_peaks.bed` - test peaks

See the README.Rmd to follow steps that will reproduce
the figure.

