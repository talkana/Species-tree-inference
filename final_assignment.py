#!/usr/bin/env python3
# Natalia Rutecka - Comparative genomics final assignment
# Code for downloading genomes and genes is based on code from Coimbra et al 2020 (https://github.com/UdeS-CoBIUS/EXECT)

import argparse
from Bio import SeqIO, AlignIO, Entrez, GenBank
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor, ParsimonyTreeConstructor, \
    ParsimonyScorer, NNITreeSearcher
import subprocess
from Bio.Phylo.Consensus import bootstrap
import multiprocessing
import os
import shutil


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Builds consensus tree and supertree for a collection of species' genome names.")
    parser.add_argument("-n", required=True,
                        help="path to a file containing species' genome names")
    parser.add_argument('--bootstrap', dest="b", action='store_true', default=False,
                        help="Calculate bootstrap supports for gene trees")
    parser.add_argument('-f', default=0, type=float, required=False,
                        help='Mean bootstrap support threshold (use with --bootstrap). Will not use gene trees with mean bootstrap support < f')
    parser.add_argument('-t', default=1, type=int, required=False,
                        help='Number of threads to use')
    parser.add_argument('-o', required=True,
                        help='Name of the output folder')
    parser.add_argument('-oc', required=False, help="Name of file containing gene clusters",
                        default="clustered_sequences")
    parser.add_argument('-ot', required=False, help="Name of file containing gene trees",
                        default="gene_trees")
    parser.add_argument('-os', required=False, help="Name of file containing all gene sequences",
                        default="all_sequences")
    parser.add_argument('-ct', required=False, help="Majority consensus threshold",
                        default=0.2)
    parser.add_argument('-e', required=False, help="Your e-mail address (for Entrez download)",
                        default="n.rutecka@student.uw.edu.pl")
    parsed = parser.parse_args()
    return parsed


def process_genome_list(file):
    id_list = []
    new_names = open(file)
    for line in new_names.readlines():
        genome_id = line.split('_')
        genome_id = genome_id[-2] + '_' + genome_id[-1].strip('\n')
        id_list.append(genome_id)
    return id_list


def download_genomes(genome_list, folder, email):
    Entrez.email = email
    for genome_id in genome_list:
        filename = folder + "/DB_" + genome_id + ".gbk"
        handle = Entrez.efetch(db='nuccore', id=genome_id, rettype="gbwithparts", retmode="text")
        out_handle = open(filename, "w")
        for line in handle:
            out_handle.write(line)
        out_handle.close()


def download_genes(genome_path, genes_path):
    cds_file = open(genes_path, "w")
    for file in os.listdir(genome_path):
        parser = GenBank.RecordParser()
        rc_file = open(genome_path + "/" + file)
        record = parser.parse(rc_file)
        rc_file.close()
        for feature in record.features:
            if feature.key == 'CDS':
                pid = ""
                trans = ""
                for qualifier in feature.qualifiers:
                    if qualifier.key == '/locus_tag=':
                        pid = qualifier.value[1:-1]
                        pid = pid.replace('"', '')
                    if qualifier.key == '/translation=':
                        trans = qualifier.value[1:-1]
                cds_file.write(">%s|%s" % (record.locus, pid + "\n"))
                cds_file.write("%s" % trans + "\n")
    cds_file.close()
    shutil.rmtree(genome_path)


def cluster_genes(seq_filename, results_filename, threads):
    subprocess.call(
        f"OMP_NUM_THREADS={threads} mmseqs easy-cluster {seq_filename} {results_filename} {results_filename}_tmp -c 0.5",
        shell=True)


def is_bijective(tree, taxa):
    for taxon in taxa:
        if tree.count(taxon) != 1:
            return False
    return True


def get_mmseq_clusters(mmseq_folder):
    file = open(mmseq_folder + "_cluster.tsv").read().split("\n")
    clusters = []
    for i in range(0, len(file) - 1):
        ids = file[i].split()
        for cluster in clusters:
            if ids[0] in cluster or ids[1] in cluster:
                cluster.add(ids[0])
                cluster.add(ids[1])
                break
        else:
            new_cluster = set()
            new_cluster.add(ids[0])
            new_cluster.add(ids[1])
            clusters.append(new_cluster)

    informative = []
    for cluster in clusters:
        if len(cluster) > 3:
            informative.append(cluster)
    print(f"Mamy {len(informative)} informatywnych klastrów")

    return informative


def align_genes(gene_file, aln_file):
    subprocess.call(f"mafft --anysymbol {gene_file} > {aln_file}", shell=True)


def has_positive_lengths(tree):
    i = 0
    newick = tree.format(fmt="newick")
    while i < len(newick):
        if newick[i] == ":":
            i += 1
            st = i
            while newick[i].isalnum() or newick[i] in ".-+":
                i += 1
            bl = float(newick[st:i])
            if bl < 0:
                return False
        else:
            i += 1
    return True


def rename_to_taxa(newick):
    i = 0
    cleaned = ""
    newick = newick.replace("-", "")
    while i < len(newick):
        if newick[i].isalnum() or newick[i] in "/_":
            if newick[i:i + 5] != "Inner":
                j = i
                while j < len(newick) and (newick[j].isalnum() or newick[j] in "/_"):
                    j += 1
                cleaned += newick[i:j]
                i = j
            while newick[i].isalnum() or newick[i] in "/_|":
                i += 1
        else:
            cleaned += newick[i]
            i += 1
    return cleaned


def nj_tree(aln):
    calculator = DistanceCalculator('blosum62')
    distance_matrix = calculator.get_distance(aln)
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(distance_matrix)
    return tree


def parsimony_tree(aln):
    scorer = ParsimonyScorer()
    searcher = NNITreeSearcher(scorer)
    constructor = ParsimonyTreeConstructor(searcher)
    tree = constructor.build_tree(aln)
    return tree


def get_supports(aln, btrees_path, method, bn=50):
    aln_reps = bootstrap(aln, bn)
    btrees_file = open(btrees_path, "w")

    for aln in aln_reps:  # save bootstrap trees
        if method == "nj":
            tree = nj_tree(aln)
        elif method == "parsimony":
            tree = parsimony_tree(aln)
        else:
            raise ValueError(f"Inference method {method} is not supported")
        nwk = tree.format(fmt="newick")
        btrees_file.write(nwk)
    btrees_file.close()


def infer_tree(aln_path, bstr, btrees_path, rename=True):
    method = "nj"
    aln = AlignIO.read(aln_path, 'fasta')
    tree = nj_tree(aln)

    if not has_positive_lengths(tree):
        method = "parsimony"
        tree = parsimony_tree(aln)

    if bstr:
        get_supports(aln, btrees_path, method)
    nwk = tree.format(fmt="newick")
    if rename:
        nwk = rename_to_taxa(nwk)
    return nwk


def get_maj_consensus(thres, trees_path, consensus_path, threads):
    subprocess.call(
        f"sumtrees.py --force-unrooted -o {consensus_path} -F newick --suppress-annotations -m {threads} -s consensus -f {thres} -d2 {trees_path}",
        shell=True)


def get_supertree(gene_trees_path):
    subprocess.call(["fasturec", "-G", gene_trees_path, "-Y"])


def names_to_seqs(seqfile):
    seqs = {}  # map from name to seq
    file = open(seqfile)
    sequences = SeqIO.parse(file, "fasta")
    for seq in sequences:
        seqs[seq.id] = seq.seq
    file.close()
    return seqs


def get_mean_support(newick):
    s = 0
    nr = 0
    i = 0
    while i < len(newick):
        if newick[i] == ")":
            if newick[i + 1] == ";":
                i += 1
                continue
            st = i + 1
            while i < len(newick) and newick[i] != ":":
                i += 1
            bl = float(newick[st:i])
            s += bl
            nr += 1
        else:
            i += 1
    perc = s / nr
    return perc


def save_gene_tree(tree, trees_path, mean_bootstrap):
    trees_file = open(trees_path, "a")
    if mean_bootstrap and get_mean_support(tree) >= mean_bootstrap:
        filtered_path = f"{trees_path}_filtered_{mean_bootstrap}"
        filtered_file = open(filtered_path, "a")
        filtered_file.write(f"{tree}")
        filtered_file.close()
    trees_file.write(tree)
    trees_file.close()


def save_sequences(sequences_name, clusters_name, outfolder, seqname):
    names = []
    clusters = get_mmseq_clusters(clusters_name)
    name_to_seq = names_to_seqs(sequences_name)
    for i in range(len(clusters)):
        subfolder = f"{outfolder}/{i + 1}"
        os.mkdir(subfolder)
        seq_name = f"{subfolder}/{seqname}"
        names.append(seq_name)
        seqfile = open(seq_name, "w")
        for seq_id in clusters[i]:
            seq = name_to_seq[seq_id]
            new_seq = f'>{seq_id}\n{seq}\n'
            seqfile.write(new_seq)
        seqfile.close()
    return names


def run_single_inference(seq, nr, trees_name, bstr, mean_support):
    print(f"Przetwarzam klaster {nr}")
    aln_name = seq + "_aligned"
    gtree_name = seq + "_gtree"
    gtree_supports_name = seq + "_gtree_supports"
    btrees_path = seq + "_btrees"
    align_genes(seq, aln_name)
    if bstr:
        gtree = infer_tree(aln_name, bstr, btrees_path, False)
    else:
        gtree = infer_tree(aln_name, bstr, btrees_path)
    gfile = open(gtree_name, "w")
    gfile.write(f"{gtree};")
    gfile.close()
    if bstr:
        subprocess.call(
            f"sumtrees.py -d2 -o {gtree_supports_name} -F newick --suppress-annotations -t {gtree_name} {btrees_path}",
            shell=True)
        gtree = open(gtree_supports_name).read().replace("[&U] ", "")
        gtree = rename_to_taxa(gtree)
    save_gene_tree(gtree, trees_name, mean_support)


def run_inference(seqs, trees_name, bst, mean_support, threads):
    seqs_labelled = [(seqs[i], i + 1, trees_name, bst, mean_support) for i in range(len(seqs))]
    pool = multiprocessing.Pool(threads)
    pool.starmap_async(run_single_inference, seqs_labelled)
    pool.close()
    pool.join()


def save_bijective(taxa, tree_file, bij_tree_file):
    bij_file = open(bij_tree_file, "w")
    with open(tree_file) as t:
        trees = t.read().split()
        for tree in trees:
            if is_bijective(tree, taxa):
                bij_file.write(f"{tree}\n")
        bij_file.close()


def main():
    # set all filenames
    args = parse_arguments()
    os.mkdir(args.o)
    sequences_name = args.os
    clusters_name = f"{args.o}/{args.oc}"
    trees_name = f"{args.o}/{args.ot}"
    consensus_path = f"{args.o}/{args.ot}"
    seqname = "seq"
    bij_trees_name = trees_name + "_bijective"
    genomes_path = args.o + "/genomes"  # temporary
    os.mkdir(genomes_path)

    # run pipeline
    idlist = process_genome_list(args.n)
    download_genomes(idlist, genomes_path, args.e)
    download_genes(genomes_path, sequences_name)
    cluster_genes(sequences_name, clusters_name, args.t)
    seqs = save_sequences(sequences_name, clusters_name, args.o, seqname)
    run_inference(seqs, trees_name, args.b, args.f, args.t)
    save_bijective(idlist, trees_name, bij_trees_name)
    get_supertree(trees_name)
    get_maj_consensus(args.ct, bij_trees_name, consensus_path, args.t)


if __name__ == "__main__":
    main()
