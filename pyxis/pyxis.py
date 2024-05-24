#!/usr/bin/env python

"""
Command-line script for finding enriched motifs in genomic peak regions

Analogous to HOMER's findMotifsGenome.pl command
"""
from . import myutils as myutils
from pyxis import __version__
import argparse
import numpy as np
import os
import pyfaidx
import pandas as pd
import seqlogo
import sys


def main():
    parser = argparse.ArgumentParser(
        prog="pyxis",
        description="Pyxis, a command-line tool for finding enriched motifs in genomic peak regions",
        usage='%(prog)s peaks.bed refgenome.fa motif.pwms [options]'
    )

    # ---------------------------------- Define and parse arguments --------------------------------------
    # Positional arguments (required)
    parser.add_argument("peaks_bed", help="BED file of ChIP peaks", metavar="PEAKS", type=str)
    parser.add_argument("ref_fasta", help="FASTA reference genome file", metavar="REF", type=str)
    parser.add_argument("motifs_pwms", help="PWMS file specifying PWMs of known motifs", 
                        metavar="PWMS", type=str)

    # Optional arguments
    parser.add_argument("-o", "--out", help="Write output to directory. Default: stdout",
                        metavar="DIR", type=str, required=False)
    parser.add_argument("-b", "--background", help="Specify own background nucleotide frequencies."
                        " Default: randomly generated from reference genome", metavar="BKGD", type=str, required=False)
    parser.add_argument("-p", "--pval", help="Use specified p-value for enrichment significance. Default: 1e-5", 
                        metavar="VAL", type=int, required=False)
    parser.add_argument("-s", "--seqlogo", help="Generate sequence logo for enriched motif", 
                        action="store_true", required=False)
    parser.add_argument("--version", help="Print version and exit", action="version",
                        version='{version}'.format(version=__version__))

    # Parse arguments
    args = parser.parse_args()

    # ------------------------- Set up output file -----------------------------------
    if args.out is None:
        outf = sys.stdout
    else: outf = open(args.out, "w")

    # --------------------------- Loading in input files -------------------------------
    # load FASTA file (indexed with pyfaidx)
    if args.ref_fasta is not None:
        if not os.path.exists(args.ref_fasta):
            myutils.ERROR("{fasta} does not exist".format(fasta=args.ref_fasta))
        reffasta = pyfaidx.Fasta(args.ref_fasta)
    else:
        myutils.ERROR("FASTA file was not specified.")

    # load foreground peaks BED file
    if args.peaks_bed is not None:
        if not os.path.exists(args.peaks_bed):
           myutils.ERROR("{foreground} does not exist".format(foreground=args.peaks_bed))
        f_seq, total_f_peaks = myutils.ReadBED(args.peaks_bed, reffasta)
    else: 
        myutils.ERROR("Foreground peaks BED file not specified.")

    # load pwms file
    if args.motifs_pwms is not None:
        if not os.path.exists(args.motifs_pwms):
            myutils.ERROR("{motifs} does not exist".format(motifs=args.motifs_pwms))
        PWMList, pwm_names = myutils.ReadPWMS(args.motifs_pwms)
    else:
        myutils.ERROR("PWMS file was not specified.")

    # load background peaks BED file and compute nucleotide frequency for background
    if args.background is not None:
        if not os.path.exists(args.background):
           myutils.ERROR("{background} does not exist".format(background=args.background))
        b_seq, total_b_peaks = myutils.ReadBED(args.background, reffasta)
        #b_freqs = myutils.ComputeNucFreqs(b_seq)
    else: 
        outf.write("Background not specified. Generating background sequences from reference genome...")
        #b_freqs = myutils.ComputeNucFreqs(myutils.RandomBkSequence(f_seq_lens, len(f_seq), reffasta)[0])
        f_seq_lens = []
        for i in range(len(f_seq)):
            f_seq_lens.append(len(f_seq[i]))
        b_seq, total_b_peaks = myutils.RandomBkSequence(f_seq_lens, len(f_seq), reffasta)

    # ------------------- Determing score thresholds for sequences to match PWM ----------------------------
    numsim = 10000
    null_pval = 0.01
    if args.pval is not None:
        thresh_pval = args.pval
    else:
        thresh_pval = 1e-5
    
    outf.write("\nCalculating thresholds for each PWM...")
    pwm_thresholds = [] # score thresholds for each PWM
    freqs = myutils.ComputeNucFreqs(f_seq + b_seq + [myutils.ReverseComplement(item) for item in f_seq] \
                                + [myutils.ReverseComplement(item) for item in b_seq])
    for i in range(len(PWMList)):
        null_scores = [myutils.ScoreSeq(PWMList[i], 
                                myutils.RandomSequence(PWMList[i].shape[1], freqs)) for j in range(numsim)]
        thresh = myutils.GetThreshold(null_scores, null_pval)
        pwm_thresholds.append(thresh)
        outf.write("\n[" + str((i + 1)) + "/" + str(len(PWMList)) + "] .... ")
    outf.write("Done.\n")

    # -------------------------------- Test for motif enrichment ------------------------------------
    outf.write("\nPerforming motif enrichment...")
    pvals = []
    num_peak_passes = []
    num_bg_passes = []
    for i in range(len(PWMList)):
        pwm = PWMList[i]
        thresh = pwm_thresholds[i]
        num_peak_pass = np.sum([int(myutils.FindMaxScore(pwm, seq)>thresh) for seq in f_seq])
        num_peak_passes.append(num_peak_pass)
        num_bg_pass = np.sum([int(myutils.FindMaxScore(pwm, seq)>thresh) for seq in b_seq])
        num_bg_passes.append(num_bg_pass)
        pval = myutils.ComputeEnrichment(total_f_peaks, num_peak_pass, total_b_peaks, num_bg_pass)
        pvals.append(pval)
        outf.write("\n[" + str((i + 1)) + "/" + str(len(PWMList)) + "] .... ")
    outf.write("Done.\n")

    # ------------------------------ Set up output file(s) ------------------------------------------
    outf.write("\nCreating output file 'enrichment_results.tsv'...")
    motif_names = np.array(pwm_names)
    pvals_arr = np.array(pvals)
    log_pvals = np.round(np.log10(pvals), 5)
    num_peak_motifs = np.array(num_peak_passes)
    pct_peak_motifs = np.array([num_peak/total_f_peaks for num_peak in num_peak_motifs])
    num_bg_motifs = np.array(num_bg_passes)
    pct_bg_motifs = np.array([num_bg/total_b_peaks for num_bg in num_bg_motifs])
    enriched_statuses = []
    for i in range(len(PWMList)):
        outf.write("\n....")
        if (pvals[i] <= thresh_pval):
            enriched_statuses.append("Sig")
        else:
            enriched_statuses.append("Non-sig")
    results = {"motif_name": motif_names, "pval": pvals, "log_pval": log_pvals, "num_peak_motif": num_peak_motifs,
                "pct_peak_motif": pct_peak_motifs, "num_bg_motif": num_bg_motifs, "pct_bg_motif": pct_bg_motifs,
                "enriched_status": enriched_statuses}
    # put a check here later for all arrays being the same size
    output = pd.DataFrame(results)
    output = output.sort_values(by=['pval'], ascending=False)
    output.to_csv("pyxis_enrichments.tsv", sep="\t")
    outf.write(" Done.")

    # generating seqlogos for motifs
    if args.seqlogo:
        outf.write("\nCreating seqlogos...")
        for i in range(len(PWMList)):
            seq_pwm = seqlogo.Pwm(PWMList[i].transpose())
            print(seq_pwm)
            #Convert to ppm needed for plotting
            seq_ppm = seqlogo.Ppm(seqlogo.pwm2ppm(seq_pwm))
            seqlogo.seqlogo(seq_ppm, ic_scale=True, format='png', size='medium', filename=pwm_names[i]+'/'+pwm_names[i]+'.png')
            outf.write("[" + (i + 1) + "/" + len(PWMList) + "]")
        outf.write("\nDone.")
    outf.write("\n\nPyxis has run successfully, please check 'pyxis_enrichments.tsv' for motif enrichment info!\n\n")

    reffasta.close()
    outf.close()
    sys.exit(0)

if __name__ == "__main__":
    main()
