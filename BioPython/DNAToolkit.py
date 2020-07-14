import collections
from structures import *


def validateSequence(dna):
    tmpseq = dna.upper()
    for nuc in tmpseq:
        if nuc not in Nucleotides:
            return False
    return tmpseq

def countNucFrequency(dna):
    tmpFreqDict = {"A": 0, "C": 0, "G": 0, "T": 0}
    for nuc in dna:
        tmpFreqDict[nuc] += 1
    return tmpFreqDict

def countNucFrequencyCollections(dna):
    return dict(collections.Counter(dna))

def transcription(dna):
    #DNA -> RNA Transcription
    """DNA -> RNA Transcription"""
    return dna.replace("T", "U")

def reverse_complement(dna):
    return ''.join([DNA_ReverseComplement[nuc] for nuc in dna])[::-1]
    #Python way of doing
    #mapping = str.maketrans('ATCG', 'TAGC')
    #return seq.translate(mapping)[::-1]

def gc_content(dna):
    """GC Content in a DNA/RNA sequence"""
    return round((dna.count('C') + dna.count('G')) / len(dna) * 100)

def gc_content_subsec(dna, k=20):
    """GC Content in a DNA/RNA sub-sequence length k. k=20 by default"""
    res = []
    for i in range(0, len(dna) - k + 1, k):
        subseq = dna[i:i + k]
        res.append(gc_content(subseq))
    return res

def translate_seq(dna, init_pos = 0):
    """Translates a DNA sequence into an aminoacid sequence"""
    return [DNA_Codons[dna[pos:pos+3]] for pos in range(init_pos, len(dna) - 2, 3)]

def codon_usage(dna, aminoacid):
    tmpList = []
    for i in range(0, len(dna) - 2, 3):
        if DNA_Codons[dna[i:i+3]] == aminoacid:
            tmpList.append(dna[i:i+3])

    freqDict = dict(collections.Counter(tmpList))
    totalWight = sum(freqDict.values())
    for dna in freqDict:
        freqDict[dna] = round(freqDict[dna] / totalWight, 2)
    return freqDict

def gen_reading_frames(dna):
    frames = []
    frames.append(translate_seq(dna, 0))
    frames.append(translate_seq(dna, 1))
    frames.append(translate_seq(dna, 2))
    
    frames.append(translate_seq(reverse_complement(dna), 0))
    frames.append(translate_seq(reverse_complement(dna), 1))
    frames.append(translate_seq(reverse_complement(dna), 2))
    return frames

def proteins_from_rf(amino_acid_seq):
    """Compute all possible proteins in an aminoacid seq anf return a list of possible proteins"""
    current_prot = []
    proteins = []
    for amino_acid in amino_acid_seq:
        if amino_acid == "_":
            # STOP accumulating amino acids if _ - STOP was found
            if current_prot:
                for p in current_prot:
                    proteins.append(p)
                current_prot = []
        else:
            # START accumulating amino acids if M - START was found
            if amino_acid == "M":
                current_prot.append("")
            for i in range(len(current_prot)):
                current_prot[i] += amino_acid
    return proteins

def all_proteins_from_open_reading_frames(seq, startReadPos=0, endReadPos=0, ordered=False):
    """Compute all possible proteins from all open reading frames"""
    """Protein Serch DB?"""
    """API for pulling protein info"""
    if endReadPos > startReadPos:
        rfs = gen_reading_frames(seq[startReadPos: endReadPos])
    else:
        rfs = gen_reading_frames(seq)

    res = []
    for rf in rfs:
        prots = proteins_from_rf(rf)
        for p in prots:
            res.append(p)
    
    if ordered:
        return sorted(res, key=len, reverse=True)
    return res

