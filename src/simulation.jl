#generate amino acid alignments
function simAlignments(numFiles,numSeqs,parentdir)
	infiles = Vector{String}(undef, numFiles)
	for f in 1:numFiles
		file, io = mktemp(parentdir)
    lenSeq = rand(100:1000)
		infiles[f] = file
    w = FASTX.FASTA.Writer(io)
		for n in 1:numSeqs
			seq = BioSequences.randaaseq(lenSeq) # all seqs need to have same length
			rec = FASTX.FASTA.Record(string(n),seq)
			write(w, rec)
		end
    close(io)
	end
	return(infiles)
end