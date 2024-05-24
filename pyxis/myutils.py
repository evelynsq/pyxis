"""
Utilities for pyxis
"""
import pyfaidx
import random
import numpy as np
import pandas as pd
import scipy.stats
import sys

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def ERROR(msg):
	"""
	Print an error message and die

	Parameters
	----------
	msg : str
	   Error message to print
	"""
	sys.stderr.write(bcolors.FAIL + "[ERROR]: " + bcolors.ENDC + "{msg}\n".format(msg=msg) )
	sys.exit(1)

# ----------------------------- Read in input files -----------------------------------------
def ReadBED(peaks, refgenome):
    """ 
    Load in + get info on peaks BED file

    Parameters
    ----------
    peaks : file 
        BED of genomic peak locations
    refgenome : pyfaidx object
        indexed reference genome fasta

    Returns
    -------
    seq : list of peak sequences
    total_peaks : total number of peak sequences
    """
    seq = [] 
    with open(peaks, 'r') as bed:
        for line in bed:
            curr = line.strip().split("\t")
            if (len(curr) != 6):
                ERROR("BED file incorrectly formatted, please verify format matches what is specified in README.md.")
            if (curr[0] != "chr"):
                seq.append(WriteFastaSeq(refgenome, "chr"+curr[0], int(curr[1]), int(curr[2])))
    total_peaks = len(seq)
    return seq, total_peaks

def ReadPWMS(motifs):
    """ 
    Load in + get info on peaks BED file

    Parameters
    ----------
    peaks : file 
        BED of genomic peak locations

    Returns
    -------
    PWMList : list of peak sequences
    pwm_names : total number of peak sequences
    """
   # pwmfile = pd.read_csv(motifs, sep='\t')
    PWMList = []
    pwm_names = []
    with open(motifs, 'r') as pwms:
        currPWM = []
        for line in pwms:
            currline = line.strip().split("\t")
            if (len(currline) == 3):
                pwm_names.append(currline[0])
            elif (len(currline) == 4):
                currPWM.append(currline)
            elif (currline == ['']):
                currPWM = (np.array(currPWM)).astype(np.float)
                PWMList.append(currPWM.transpose())
                currPWM = []
            else:
                ERROR("PWMS file incorrectly formatted, please verify format matches what is specified in README.md.")
    return PWMList, pwm_names
    
def WriteFastaSeq(refgenome, chrom, start, end):
    """
    Write out reference genome sequence for a specified region within a chromosome

    Parameters
    ----------
    refgnome : pyfaidx object
        indexed reference genome fasta
    chrom : str
        chromosome of interest
    start : int
        starting position of sequence
    end : int
        ending position of sequence

    Returns
    -------
    fastaseq : string of specified reference genome sequence
    """
    fastaseq = refgenome[chrom][start:end].seq
    return fastaseq

# ---------- Generate background sequences for comparison with peak sequences ---------------
def ComputeNucFreqs(sequences):
    """ Compute nucleotide frequencies of a list of sequences
    
    Parameters
    ----------
    sequences : list of str
       List of sequences
       
    Returns
    -------
    freqs : list of float
       Frequencies of A, C, G, T in the sequences
    """
    nucs = {"A": 0, "C": 0, "G": 0, "T": 0} # hold number of nucleotides seen in sequence 
    freqs = [0.25, 0.25, 0.25, 0.25]        # compute frequency of order A, C, G, T
    total_nucs = 0
    for i in range(len(sequences)):
        for j in range(len(sequences[i])):
            if (sequences[i][j] == 'A'):
                nucs['A'] += 1
            elif (sequences[i][j] == 'C'):
                nucs['C'] += 1
            elif (sequences[i][j] == 'G'):
                nucs['G'] += 1
            elif (sequences[i][j] == 'T'):
                nucs['T'] += 1
            total_nucs += 1
    freqs[0] = nucs['A'] / total_nucs
    freqs[1] = nucs['C'] / total_nucs
    freqs[2] = nucs['G'] / total_nucs
    freqs[3] = nucs['T'] / total_nucs
    return freqs

def RandomBkSequence(seq_lens, num_seqs, refgenome):
    """ 
    Generate random background sequences from reference genome 
    and having same length and number as given peak sequences
    
    Parameters
    ----------
    seq_lens : list
        list of sequence lengths for background sequences
    num_seqs : int
        number of background sequences to generate
    refgenome : pyfaidx object
        indexed reference genome fasta
    
    Returns
    -------
    seqs : list
       list of random background sequences from reference genome
    num_b_seqs : int
        number of background sequences
    """
    seqs = []
    chromosomes = list(refgenome.keys())
    for i in range(num_seqs):
        chr = random.choice(chromosomes)
        start = random.randint(0, len(refgenome[chr]) - seq_lens[i])
        end = start + seq_lens[i] - 1
        seqs.append(WriteFastaSeq(refgenome, chr, start, end))
    num_b_seqs = len(seqs)
    return seqs, num_b_seqs

# ------------ Obtain max PWM score for sequences and their reverse complements ---------------
def ScoreSeq(pwm, sequence):
    """ Score a sequence using a PWM
    
    Parameters
    ----------
    pwm : 2d np.array
       Position weight matrix
    sequence : str
       Sequence of nucleotides to be scored
       
    Returns
    -------
    score : float
       PWM score of the sequence
    """
    score = 0
    # Compute PWM score by adding up scores in cells of PWM that match sequence at each position
    for j in range(len(sequence)):
        if (sequence[j] == 'A'):
            score += pwm[0][j]
        elif (sequence[j] == 'C'):
            score += pwm[1][j]
        elif (sequence[j] == 'G'):
            score += pwm[2][j]
        elif (sequence[j] == 'T'):
            score += pwm[3][j]
    return score

def ReverseComplement(sequence):
    """ Get the reverse complement of a sequence
    
    Parameters
    ----------
    sequence : str
      Sequence of nucleotides
      
    Returns
    -------
    revcomp : str
      Reverse complement of sequence
    """
    revcomp = ""
    for i in range(len(sequence)):
        idx = len(sequence) - i - 1
        if (sequence[idx] == 'A'):
            revcomp += 'T'
        elif (sequence[idx] == 'C'):
            revcomp += 'G'
        elif (sequence[idx] == 'G'):
            revcomp += 'C'
        elif (sequence[idx] == 'T'):
            revcomp += 'A'
    return revcomp

def FindMaxScore(pwm, sequence):
    """ Get highest PWM match for a sequence
    
    Scan a sequence with a pwm
    Compute the highest pwm score for a given sequence
    Be sure to check the forward and reverse strands!
    
    Parameters
    ----------
    pwm : 2d np.array
       PWM matrix
    sequence : str
       Sequence of nucleotides
       
    Returns
    -------
    max_score : float
       Score of top match to the PWM
    """
    if (len(pwm) == 0):
        ERROR("PWM matrix was not found, verify matrices in PWMS file match format description in README.md.")
    if (len(sequence) == 0):
        ERROR("Sequence to compute max score for was not found. Exiting")
    max_score = -1*np.inf
    n = pwm.shape[1]
    if (len(sequence) < n):
        ERROR("Sequence length is shorter than PWM length. Exiting")
    # lists of scores. scores[i] should give the score of the substring sequence[i:i+n]
    forward_scores = [0]*(len(sequence) - n + 1)
    reverse_scores = [0]*(len(sequence) - n + 1)
    rev_seq = ReverseComplement(sequence)
    for i in range(len(sequence) - n + 1): # of possible windows for n-length substrings
        forward_scores[i] = ScoreSeq(pwm, sequence[i:i+n])
        reverse_scores[i] = ScoreSeq(pwm, rev_seq[i:i+n])
    max_score = max(max(forward_scores), max(reverse_scores))
    return max_score

def RandomSequence(n, freqs):
    """ 
    Generate a random string of nucleotides of length n
    with the given nucleotide frequences
    
    Parameters
    ----------
    n : int
       Length of random string to generate
    freqs : list of float
       List of frequencies of A, C, G, T
       
    Returns
    -------
    seq : str
       random sequence of length n with the specified allele frequencies
    """
    seq = "A"*n
    nucs = ['A', 'C', 'G', 'T']
    temp_seq = []
    for i in range(len(freqs)):
        temp_seq += (nucs[i] * round(freqs[i] * n)) # add the frequency * n number of each respective nucleotide
    random.shuffle(temp_seq)
    seq = ''.join(temp_seq[:n]) # only take the first 'n' nucleotides in the shuffled sequence
    return seq

# ------------ Compute threshold for determining if peak has match to known motif --------------
def GetThreshold(null_dist, pval):
    """ 
    Find the threshold to achieve a desired p-value
    
    Given a null distribution (list of values),
    find the threshold to achieve a desired p-value
    
    Parameters
    ----------
    null_dist : list of float
       Null distribution of scores
    pval : float
       % of null_dist that should be above the threshold returned
       
    Returns
    -------
    thresh : float
       Threshold to achieve the desired p-value    
    """
    thresh = 0 
    # score threshold to obtain p-value < pval
    thresh = np.percentile(null_dist, 100 - (pval * 100))
    return thresh

# ------------------ Test if the motifs are enriched in peaks -------------------------------------
def ComputeEnrichment(peak_total, peak_motif, bg_total, bg_motif):
    """ 
    Compute fisher exact test to test whether motif enriched in bound sequences
    
    Parameters
    ----------
    peak_total : int
       Number of total peaks
    peak_motif : int
       Number of peaks matching the motif
    bg_total : int
       Number of background sequences (peak or not)
    bg_motif : int
       Number of background sequences matching the motif
       
    Returns
    -------
    pval : float
       Fisher Exact Test p-value    
    """
    pval = -1
    # contingency table
    table = [[peak_motif, bg_motif], [peak_total - peak_motif, bg_total - bg_motif]]
    odds, pval = scipy.stats.fisher_exact(table)
    return pval
