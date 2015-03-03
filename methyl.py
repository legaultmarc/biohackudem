#!/usr/bin/env python

from __future__ import print_function
import multiprocessing
from subprocess import Popen, PIPE
import argparse
import os

from Bio import SeqIO, pairwise2
import numpy as np
from gepyto.structures.sequences import Sequence
import matplotlib.pyplot as plt
import pandas as pd


def main():
    options = parse_arguments()

    root = options.out
    fastqs = (options.fwd1, options.fwd2, options.rev1, options.rev2)

    # Step 1: Compute sequencing metrics.
    sequence_lengths = []
    for fastq in fastqs:
        with open(fastq) as f:
            for record in SeqIO.parse(f, "fastq"):
                sequence_lengths.append(len(record.seq))

    sequence_lengths = np.array(sequence_lengths)
    mean_seq_length = int(round(np.mean(sequence_lengths)))
    std_seq_length = np.std(sequence_lengths)
    

    # Step 2: Merge the reads for reverse and forward.
    # This step is dependant on the FLASH tool.
    step2_dir = os.path.join(root, "merged")
    if not os.path.isdir(step2_dir):
        os.makedirs(step2_dir)

    pairs = ((options.fwd1, options.fwd2), (options.rev1, options.rev2))
    for pe1, pe2 in pairs:
        mode = "forward" if pe1 == options.fwd1 else "reverse"

        proc = Popen([
            options.flash, pe1, pe2,
            "--read-len={}".format(mean_seq_length),
            "--output-dir={}".format(step2_dir),
            "--output-prefix={}".format(mode),
            "--threads={}".format(multiprocessing.cpu_count() - 1)
        ], stdout=PIPE, stderr=PIPE)

        stdout, stderr = proc.communicate()

        with open(os.path.join(step2_dir, mode + "_stdout.log"), "wb") as f:
            f.write(stdout)
        with open(os.path.join(step2_dir, mode + "_stderr.log"), "wb") as f:
            f.write(stderr)

    flash_fwd = os.path.join(step2_dir, "forward.extendedFrags.fastq")
    flash_rev = os.path.join(step2_dir, "reverse.extendedFrags.fastq")

    assert os.path.isfile(flash_fwd), "File not found. Check FLASH output."
    assert os.path.isfile(flash_rev), "File not found. Check FLASH output."

    # Step 3: Velvet
    # Run the velveth algorithm.
    step3_dir = os.path.join(root, "assembly")
    if not os.path.isdir(step3_dir):
        os.makedirs(step3_dir)

    for mode, fastq in (("fwd", flash_fwd), ("rev", flash_rev)):
        cur_dir = os.path.join(step3_dir, mode)
        if not os.path.isdir(cur_dir):
            os.makedirs(cur_dir)

        proc = Popen([
            options.velveth, cur_dir, str(options.velvet_k),
            "-fastq", fastq,
        ], stdout=PIPE, stderr=PIPE)

        stdout, stderr = proc.communicate()

        log = os.path.join(cur_dir, "velveth_stdout.log")
        with open(log, "wb") as f:
            f.write(stdout)
        log = os.path.join(cur_dir, "velveth_stderr.log")
        with open(log, "wb") as f:
            f.write(stderr)

        # Run the velvetg tool.
        proc = Popen([
            options.velvetg, cur_dir
        ], stdout=PIPE, stderr=PIPE)

        stdout, stderr = proc.communicate()

        with open(os.path.join(cur_dir, "velvetg_stdout.log"), "wb") as f:
            f.write(stdout)
        with open(os.path.join(cur_dir, "velvetg_stderr.log"), "wb") as f:
            f.write(stderr)

    contigs_fwd = os.path.join(step3_dir, "fwd", "contigs.fa")
    contigs_rev = os.path.join(step3_dir, "rev", "contigs.fa")

    # Step 4: Compare contigs.
    sequences, dists = compare_contigs(contigs_fwd, contigs_rev, options.bbc_k,
                                       options.plot_contig_distances)


    dists = dists.sort("dist")
    cutoff = dists["dist"].mean() - dists["dist"].std()
    paired_contigs = []
    for row in dists.iterrows():
        if row.dist > cutoff:
            break

        paired_contigs.append((row.contig1, row.contig2))
    
    methyl_found = 0
    # We need to align paired contigs to identify variation.
    for contig1, contig2 in paired_contigs:
        seq1 = sequences["fwd"][contig1]
        seq2 = sequences["rev"][contig2]
        alignments = pairwise2.align.localxx(seq1.seq, seq2.seq)
        start, end = alignments[3], alignments[4]

        for pos in xrange(start,end+1):
            if alignments[0][pos]=='C' and alignments[1][pos]=='C' or
                    alignments[1][pos]=='A' and alignments[0][pos]=='A':


def compare_contigs(fwd, rev, k, plot=False):
    """Compare fastq contigs.

    Takes a path to the forward and to the reverse fastq contigs.

    """
    bbcs = {
        "fwd": {},
        "rev": {},    
    }

    sequences = {
        "fwd": {},
        "rev": {},
    }

    for mode, fastq in (("fwd", fwd), ("rev", rev)):
        with open(fastq) as f:
            for record in SeqIO.parse(f, "fasta"):
                 # Compute the bbc for this record.
                 seq = Sequence(record.id, record.seq, "DNA")
                 sequences[mode][record.id] = seq
                 a = set(list("ATGC"))
                 bbcs[mode][record.id] = seq.bbc(k, alphabet=a).reshape(1, 16)


    distances = []
    for contig1, bbc1 in bbcs["fwd"].items():
        for contig2, bbc2 in bbcs["rev"].items():
            distances.append((contig1, contig2, np.sum((bbc1 - bbc2) ** 2)))
    distances = pd.DataFrame(distances, columns=["contig1", "contig2", "dist"])

    if plot:
        dists = distances["dist"].values
        plt.hist(dists, bins=40, edgecolor="white", facecolor="black")
        plt.show()

    return sequences, distances


def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument("--fwd1",
        help="Forward read 1.", type=str, required=True,
    )

    parser.add_argument("--fwd2",
        help="Forward read 2.", type=str, required=True,
    )

    parser.add_argument("--rev1",
        help="Reverse read 1.", type=str, required=True,
    )

    parser.add_argument("--rev2",
        help="Reverse read 2.", type=str, required=True,
    )

    parser.add_argument("--out", "-o",
        help="Work directory.", type=str, default="methylation_analysis"
    )

    parser.add_argument("--flash",
        help="Path to the flash binary.", type=str, default="flash"
    )

    parser.add_argument("--velveth",
        help="Path to the velveth binary.", type=str, default="velveth"
    )

    parser.add_argument("--velvetg",
        help="Path to the velvetg binary.", type=str, default="velvetg"
    )

    parser.add_argument("--velvet-k",
        help="The velvet hash length.", type=int, default=30
    )

    parser.add_argument("--bbc-k",
        help="k radius for the BBC.", type=int, default=20
    )

    parser.add_argument("--plot-contig-distances",
        help="Plot pairwise contig distances.", action="store_true"
    )

    return parser.parse_args()


if __name__ == "__main__":
    main()
