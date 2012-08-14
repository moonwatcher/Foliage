Foliage is a software library for efficiently constructing ancestral recombination graphs in memory, written in Java.
It facilitates the implementation of various algorithms and statistical analysis by providing a [coherent API](http://moonwatcher.github.com/Foliage/) for the graphs.
It is capable of cropping the graphs along the genome coordinates, drawing [Graphviz](http://www.graphviz.org) diagrams of the graphs 
and correlation heat maps that measure linkage across the genome by comparing the distance between the marginal trees.

Here are some graphs for abour 16 SPNs on the MHC region. Graphs can be constructed as [bi furcated graphs](http://moonwatcher.github.com/Foliage/HLA_2135_2150.pdf), which corresponds to the classic Wright-Fisher model, 
or [multi furcated graphs](http://moonwatcher.github.com/Foliage/HLA_mf_2135_2150.pdf), which can be significantly smaller, and collapse events that have no direct evidence in the dataset.

For the sake of demonstration, here are graphs spanning some 2150 SNPs on the MHC: [bi furcated](http://moonwatcher.github.com/Foliage/HLA.pdf) and [multi furcated](http://moonwatcher.github.com/Foliage/HLA_mf.pdf).
Those are quite big but do a good job of demonstrating the difference in size when allowing the multifurcated reduction.
For the purpose of mesuring pairwise correlations the multifurcated version can be argued to be more accuret as well as faster to work with.

 
The integrated [margarita](http://www.sanger.ac.uk/Software/analysis/margarita) library, written by Mark Minichiello, 
can be used to infer plausible ARGs from fastPHASE haplotype sequence files.

Foliage is released under the terms of the [GNU GPL](http://www.gnu.org/licenses/old-licenses/gpl-2.0.html).
It was written during my work period at the [Richard Durbin](http://www.sanger.ac.uk/Teams/faculty/durbin) lab 
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
