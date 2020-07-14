from DNAToolkit import *
from utilities import colored
import random

randomStr = ''.join([random.choice(Nucleotides) for nuc in range(50)])

DNAStr = validateSequence(randomStr)

testDNAStr = "AThTCGT"

print(f'\nSequence: {colored(DNAStr)}\n')
print(f'[1] + Sequence Lenght: {len(DNAStr)}\n')
print(colored(f'[2] + Nucleotide Frequency {countNucFrequency(DNAStr)}\n'))
print(colored(f'[3] + DNA -> RNA Transcription {transcription(DNAStr)}\n'))
print(colored(f"[4] + DNA String + Complement + Reverse Complement: \n5' {DNAStr} 3'"))
print(f"   {''.join(['|' for c in range(len(DNAStr))])}")
print(colored(f"3' {reverse_complement(DNAStr)} 5'\n"))
print(colored(f"5' {reverse_complement(DNAStr[::-1])} 3' [Rev Complement]"))
print(colored(f"[5] + GC Content: {gc_content(DNAStr)}%\n"))
print(f"[6] + GC Content in Subsection k=5: {gc_content_subsec(DNAStr, k=5)}\n")
print(f"[7] + Aminoacids Sequence from DNA: {translate_seq(DNAStr, 0)}\n")
print(f'[8] + Codon Frequency (L): {codon_usage(DNAStr, "L")}\n')
print(f'[9] + Reading_frames:')
for frame in gen_reading_frames(DNAStr):
    print(frame)

test_rf_frame = ['L', 'M', 'T', 'A', 'L', 'V', 'V', 'L', 'S', 'R', 'R', 'G', 'S', '_', 'G', 'H']
print('\n[10] + All proteins in 6 open reading frames:')
for prot in all_proteins_from_open_reading_frames(DNAStr, 0, 0, True):
    print(f'{prot}')
