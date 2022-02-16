## Phylogenetic pipeline
### Final assignment for comparative genomics classes
### author: Natalia Rutecka
The pipeline infers a phylogenetic tree for a collection of species' genome names. 
It downloads the genomes, clusters the genes into homology groups, infers a gene tree for each of them and then uses supertree and consensus methods to infer the final species tree.
### Required software
To successfully run the pipeline you need to install: 
- [MMSeq2](https://github.com/soedinglab/MMseqs2)
- [fasturec](http://bioputer.mimuw.edu.pl/gorecki/fasturec/) and add its executable to your PATH
- [Mafft](https://anaconda.org/bioconda/mafft)
- [dendropy](https://dendropy.org/programs/sumtrees.html
)

### Usage
`python3 final_assignment.py [-h] -n N [--bootstrap] [-f F] [-t T] -o O [-oc OC] [-ot OT] [-os OS] [-ct CT] [-e E]`

arguments:
-  -h, --help   show this help message and exit
-  -n N         path to a file containing species' genome names
-  --bootstrap  Calculate bootstrap supports for gene trees
-  -f F         Mean bootstrap support threshold (use with --bootstrap). Will not use gene trees with mean bootstrap support < f
-  -t T         Number of threads to use
-  -o O         Name of the output folder
-  -oc OC       Name of file containing gene clusters
-  -ot OT       Name of file containing gene trees
-  -os OS       Name of file containing all gene sequences
-  -ct CT       Majority consensus threshold (default: 0.2)
-  -e E         Your e-mail address (for Entrez download)


### Data and results
The pipeline was tested on Corynebacteria data, 
which is a subtaxon of Corynebacteriales analysed in Coimbra et al 2020. 
The data folder contains a genome list and a tree from Coimbra et al as well as the same files modified in a way that they only contain Corynebacteria data.

Using the pipeline I inferred 6 trees: 4 supertrees and 2 majority consensus trees.
All 6 trees are available in the results folder 
together with 4 sets of gene trees
used in the inference step.

Gene trees files marked as bijective contain only trees that for each taxon have exactly one leaf 
labelled by it.

Gene trees files marked as well supported contain only trees with mean bootstrap support > 0.7 (calculated from 50 boostrap repliacates).
