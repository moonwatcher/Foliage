Foliage is a software library for efficiently constructing ancestral recombination graphs in memory, written in Java.
It facilitates the implementation of various algorithms and statistical analysis by providing an [API](http://moonwatcher.github.io/Foliage/) for the graphs. It can crop the graphs along the genome coordinates and draw [Graphviz](http://www.graphviz.org) diagrams of the graphs.
Foliage can plot various [recombination correlation heat maps](http://moonwatcher.github.io/Foliage/A1B8AVG_A3B7_2K.png) that measure linkage across the genome by comparing the distance between the marginal trees using the Felsenstein tree distance metric. It can also plot Linkage Disequilibrium or [compare the two](http://moonwatcher.github.io/Foliage/0_TTCP_LD_A1B8.png).

Here are some graphs for a small region of the Human Major Histocompatibility Complex. Graphs can be constructed as the classic [bi furcated graphs](http://moonwatcher.github.io/Foliage/6_A3B7.3.1300_1340_2.png), corresponding to the Wright-Fisher model, 
or [multi furcated graphs](http://moonwatcher.github.io/Foliage/8_A3B7.3.1300_1340_2.mf.png), which can be significantly smaller, and collapse events that have no direct evidence in the dataset.

For the sake of demonstration, here are graphs spanning some 2150 SNPs on the MHC: [bi furcated](http://moonwatcher.github.io/Foliage/HLA.pdf) and [multi furcated](http://moonwatcher.github.io/Foliage/HLA_mf.pdf).
Those are quite big but do a good job of demonstrating the difference in size when allowing the multi-furcated reduction.
For the purpose of measuring pairwise correlations the multi-furcated version can be argued to be more accurate as well as faster to work with.

The integrated [margarita](https://doi.org/10.1086/508901) library, written by Mark Minichiello, 
can be used to infer plausible ARGs from fastPHASE haplotype sequence files.

Foliage is released under the terms of the [GNU GPL](http://www.gnu.org/licenses/old-licenses/gpl-2.0.html).
It was written during my work period at the [Richard Durbin](https://www.sanger.ac.uk/person/durbin-richard/) lab 
at the [Sanger Institute](http://www.sanger.ac.uk) between November 2006 to February 2008.

An __Ancestral Recombination Graph__ is a [directed acyclic graph](http://en.wikipedia.org/wiki/Directed_acyclic_graph)
composed of the trees that relate individuals in a population at each position on the genome.
The trees vary as one moves along the chromosome, because of ancestral recombination events.


The [Apache Solr](http://lucene.apache.org/solr/) [OpenBitSet](http://lucene.apache.org/java/2_4_0/api/org/apache/lucene/util/OpenBitSet.html) 
is used as the underlying binary vector representation class.

Margarita depends on the following libraries:

 * Colt: http://acs.lbl.gov/~hoschek/colt/
 * JSci: http://jsci.sourceforge.net/

Foliage further depends on:

 * Jakarta ORO: http://jakarta.apache.org/oro/

These should be placed in the lib directory and named according to the ant build file.
