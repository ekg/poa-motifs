# building motifs using partial order alignment

## preparing the sequences

We convert the "queries" from Manuel into a FASTA file.

```
cat merged_cluster_tomtom_results/tomtom_merged_clusts_euclidean.tsv | tail -n+2 | cut -f 8 \
    | awk '{ print length($1), $0 }' | sort -nr \
    | awk '{ print $2 }' | awk '{ print ">seq_"NR; print $1; }' \
    | head -286 >tomtom_merged_clusts_euclidean.fa
```

## graph construction with spoa

This constructs a GFA-format graph representing the partial order alignment of the input sequences (`-G`), and includes the consensus sequence in the graph (`-C`).
Each sequence is aligned to the graph in order of length, adding to it progressively.
If the reverse complement of a sequence aligns better, we include that in the graph rather than the forward complement (`-R`).

```
spoa tomtom_merged_clusts_euclidean.fa -GCR >tomtom_merged_clusts_euclidean.spoa.gfa
```

## indexing the graph with vg

We then build several index structures that `vg map` requires.

1. The xg index, which is a kind of succinct positional index of the graph.
2. The GCSA2 index, which is a generalized FM-index of the kmers of the graph
3. The GBWT, which supports haplotype matching operations used to establish a haplotype score during mapping.

```
vg convert -g -x tomtom_merged_clusts_euclidean.spoa.prune-Tm5_3x.gfa >tomtom_merged_clusts_euclidean.spoa.prune-Tm5_3x.xg
vg index \
    -g tomtom_merged_clusts_euclidean.spoa.prune-Tm5_3x.gcsa -k 5 -X 1 \
    -G tomtom_merged_clusts_euclidean.spoa.prune-Tm5_3x.gbwt -T \
    tomtom_merged_clusts_euclidean.spoa.prune-Tm5_3x.xg
```

## mapping reads with haplotype scoring

Now, we can check if the original sequences map into the graph with a good haplotype score.
Haplotype scoring is provided by an extension of the [sublinear Li and Stephens model](https://doi.org/10.1186/s13015-019-0144-9) to the graph, using the GBWT.

```
( echo sequence identity haplotype.score | tr ' ' '\t';
  vg map -x tomtom_merged_clusts_euclidean.spoa.prune-Tm5_3x.xg \
      -g tomtom_merged_clusts_euclidean.spoa.prune-Tm5_3x.gcsa \
      --gbwt-name tomtom_merged_clusts_euclidean.spoa.prune-Tm5_3x.gbwt \
      -f tomtom_merged_clusts_euclidean.fa -j \
      | jq -cr 'select(.identity != null) | [.sequence, .identity, .annotation.haplotype_score] | @tsv' 
      ) >tomtom_merged_clusts_euclidean.aligned.tsv
```

To ensure that our matching isn't random, we can simply reverse the sequences and check how well they map.
This gives us the same base content and length while disrupting the sequences.

First, reverse the sequences:\

```
cat merged_cluster_tomtom_results/tomtom_merged_clusts_euclidean.tsv | tail -n+2 | cut -f 8 | awk '{ print length($1), $0 }' | sort -nr | awk '{ print $2 }' | rev | awk '{ print ">seq_"NR; print $1; }' | head -286 >tomtom_merged_clusts_euclidean.rev.fa
```

Then align, as before:

```
( echo sequence identity haplotype.score | tr ' ' '\t';
  vg map -x tomtom_merged_clusts_euclidean.spoa.prune-Tm5_3x.xg \
      -g tomtom_merged_clusts_euclidean.spoa.prune-Tm5_3x.gcsa \
      --gbwt-name tomtom_merged_clusts_euclidean.spoa.prune-Tm5_3x.gbwt \
      -f tomtom_merged_clusts_euclidean.rev.fa -j \
      | jq -cr 'select(.identity != null) | [.sequence, .identity, .annotation.haplotype_score] | @tsv' 
      ) >tomtom_merged_clusts_euclidean.rev.aligned.tsv
```

We can then combine the results into one table for comparisons.

```
( echo group sequence identity haplotype.score; tail -n+2 tomtom_merged_clusts_euclidean.aligned.tsv | awk '{ print "fwd\t"$0; }'; tail -n+2 tomtom_merged_clusts_euclidean.rev.aligned.tsv | awk '{ print "rev\t"$0; }') | tr ' ' '\t' >tomtom_merged_clusts_euclidean.fwd+rev.aligned.tsv
```

And plot the results in R:

```R
require(tidyverse)
x <- read.delim('tomtom_merged_clusts_euclidean.fwd+rev.aligned.tsv')
ggplot(x, aes(x=identity, y=haplotype.score, color=group)) + geom_density_2d(h=c(0.1,10)) + geom_point(alpha=I(1/3))
```
