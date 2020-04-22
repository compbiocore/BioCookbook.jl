"""
  extractReads(reader::XAM.BAM.Reader, chrom::String, chromrange::UnitRange)

Extracts reads from a region as a FASTQ file. Note: This will NOT keep
track of pair information!
"""
function extractReads(reader::XAM.BAM.Reader, chrom::String,
  chromrange::UnitRange, outfile::String)
  open(FASTX.FASTQ.Writer, outfile) do w
    for record in XAM.BAM.GenomicFeatures.eachoverlap(reader, chrom, chromrange)
    	rec = FASTX.FASTQ.Record(XAM.BAM.tempname(record), XAM.BAM.sequence(record),
        XAM.BAM.quality(record))
        write(w, rec)
    end
  end
end
#TODO: outfile as gzipped file