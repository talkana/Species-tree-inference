# Script used to process files containing Corynebacteriales species and species tree
# from Coimbra et al 2020 (https://github.com/UdeS-CoBIUS/EXECT)
# Saves chosen set of species and species tree pruned to the Corynebacteria clade

import dendropy

full_tree_path = "corynebacteriales_tree_full.nwk"
full_genomes_path = "corynebacteriales_full_list.txt"
chosen_tree_path = "corynebacteria_tree.nwk"
chosen_genomes_path = "corynebacteria_list.txt"


full_tree = open(full_tree_path).read()
full_tree = dendropy.Tree.get(data=full_tree, schema="newick")
all_genomes = open("corynebacteriales_full_list.txt").read().split()
genomes_file = open(chosen_genomes_path, "w")
chosen = []
species = set()

for taxon in all_genomes:
    if "Corynebacterium" in taxon:
        sp = "".join(taxon.split("_")[:2])
        if sp not in species or "sp" in taxon:
            species.add(sp)
            genomeId = " ".join(taxon.split('_'))
            genomes_file.write(f"{taxon}\n")
            chosen.append(genomeId)
genomes_file.close()

tree_chosen = full_tree.extract_tree_with_taxa_labels(labels=chosen)
for node in tree_chosen.leaf_node_iter():
    node.taxon.label = "_".join(node.taxon.label.split(" ")[-2:])
for node in tree_chosen.postorder_node_iter():
    node.edge.length = None
    node.label = None

newick = tree_chosen.as_string("newick").replace("'", "")
chosen_tree_file = open(chosen_tree_path, "w")
chosen_tree_file.write(newick)
chosen_tree_file.close()
