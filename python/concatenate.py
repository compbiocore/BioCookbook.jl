# taken from https://gist.github.com/kgori/f0532cff6500e839cb29 and modified very slightly

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import UnknownSeq, Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
import random
import tempfile
import os
import timeit

#generate alignments
def simSeq(l):
    dna = ["A", "G", "C", "T", "-"]
    seq=''
    for i in range(l):
        seq+=random.choice(dna)
    return(seq)

def simAlignments(numFiles,numSeqs,tmpdir):
    lenSeqs = [random.randint(1,1000) for x in range(numFiles)]
    infiles = []
    for i in range(numFiles):
        file = tempfile.mkstemp(dir=tmpdir)
        os.close(file[0])
        infiles.append(file[1])
        seqs = [simSeq(lenSeqs[i]) for j in range(numSeqs)]
        msa = MultipleSeqAlignment(SeqRecord(Seq(s), id=str(j)) 
               for (j,s) in enumerate(seqs))
        AlignIO.write(msa,file[1],"fasta")
    return infiles

#def cleanup(infiles):
#    for f in infiles:
#        os.remove(f)

def concatenate(infiles,outfile):
    alignments = [AlignIO.read(open(f,"r"),"fasta") for f in infiles]

    # Get the full set of labels (i.e. sequence ids) for all the alignments
    all_labels = set(seq.id for aln in alignments for seq in aln)
    
    # Make a dictionary to store info as we go along
    # (defaultdict is convenient -- asking for a missing key gives back an empty list)
    tmp = defaultdict(list)
    
    # Assume all alignments have same alphabet
    #alphabet = alignments[0]._alphabet
    
    for aln in alignments:
        length = aln.get_alignment_length()
        
        # check if any labels are missing in the current alignment
        these_labels = set(rec.id for rec in aln)
        missing = all_labels - these_labels
        
        # if any are missing, create unknown data of the right length,
        # stuff the string representation into the tmp dict
        for label in missing:
            #new_seq = UnknownSeq(length, alphabet=alphabet)
            new_seq = UnknownSeq(length)
            tmp[label].append(str(new_seq))
        
        # else stuff the string representation into the tmp dict
        for rec in aln:
            tmp[rec.id].append(str(rec.seq))
            
    # Stitch all the substrings together using join (most efficient way),
    # and build the Biopython data structures Seq, SeqRecord and MultipleSeqAlignment
    #msa = MultipleSeqAlignment(SeqRecord(Seq(''.join(v), alphabet=alphabet), id=k) 
    msa = MultipleSeqAlignment(SeqRecord(Seq(''.join(v)), id=k)
               for (k,v) in tmp.items())

#    with open(outfile, "w") as out:
    AlignIO.write(msa,outfile,"fasta")

#tmpdir = tempfile.TemporaryDirectory()
#print(tmpdir.name)
#timeit.timeit('concatenate(infiles,outfile)',
#    setup='infiles=simAlignments(10,10,tmpdir.name),outfile=tempfile.NamedTemporaryFile(dir=tmpdir).name')

# python -m timeit -s 'import tempfile; tmpdir=tempfile.TemporaryDirectory(); from concatenate import simAlignments; infiles=simAlignments(10,10,tmpdir.name); outf=tempfile.NamedTemporaryFile().name' "from concatenate import concatenate; concatenate(infiles,outf)"
#100 loops, best of 3: 2.94 msec per loop