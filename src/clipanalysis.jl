using Gadfly
import Cairo, Fontconfig

"""
  clippedStats(reader)
  clippedStats(reader, chrom::String, chromrange::UnitRange)

Returns the number of reads that were clipped from the BAM file, or
from a region as readsClipped, readsTotal, and a DataFrame containing
all of the reads mapped to the region with the read type (soft/hard/both/not)
and mean quality of the read
"""
function clippedStats(reader::XAM.BAM.Reader, chrom::String, chromrange::UnitRange)
  type = String[]
  qual = Float64[]
  readsClipped = 0
  readsTotal = 0
  for record in XAM.BAM.GenomicFeatures.eachoverlap(reader, chrom, chromrange)
    readsTotal = readsTotal + 1
    clippedStatsHelper!(record,qual,type,readsClipped)
  end
  return readsClipped, readsTotal, DataFrames.DataFrame(type=type,qual=qual)
end

function clippedStats(reader::XAM.BAM.Reader)
  type = String[]
  qual = Float64[]
  readsClipped = 0
  readsTotal = 0
  for record in reader
    readsTotal = readsTotal + 1
    clippedStatsHelper!(record,qual,type,readsClipped)
  end
  return readsClipped, readsTotal, DataFrames.DataFrame(type=type,qual=qual)
end

function clippedStatsHelper!(record::XAM.BAM.Record, qual::Vector{Float64},type::Vector{String},
  readsClipped::Int8)
  cigar = XAM.BAM.cigar_rle(record)
  push!(qual,Statistics.mean(XAM.BAM.quality(record)))
  flag = "not"
  if findfirst(x->x==convert(BioAlignments.Operation,'S'),cigar[1]) != nothing
    flag = "soft"
  end
  if findfirst(x->x==convert(BioAlignments.Operation,'H'),cigar[1])!= nothing
    if flag == "soft"
      push!(type,"both")
    else
      push!(type,"hard")
    end
    readsClipped = readsClipped + 1
  else
    if flag == "soft"
      push!(type,"soft")
      readsClipped = readsClipped + 1
    else
      push!(type,"not")
    end
  end
  return type,qual,readsClipped
end

"""
  clippedReadQual(df::DataFrame, outfile::String="clippedReadQual.pdf")

Plots a histogram of mean read quality split by whether the read was soft-clipped,
hard-clipped, both, or not clipped and saves it to outfile.
"""
function clippedReadQual(df::DataFrames.DataFrame, outfile::String="clippedReadQual.pdf")
  p = plot(df,x="qual", color="type", Geom.histogram,
    Guide.title("Quality score distribution for clipped and not in region"),
    Guide.xlabel("Mean quality"),
    Guide.ylabel("Number of reads"))
  p |> PDF(outfile)
end

#using BenchmarkTools
#@benchmark
#reader = open(BAM.Reader, "/Users/aguang/CORE/scratch/Bliss_71991_101.dedup.rg.srtd.bam",
#  index="/Users/aguang/CORE/scratch/Bliss_71991_101.dedup.rg.srtd.bam.bai")
#df = clippedStats("Ca22chr4A_C_albicans_SC5314",726523:727338,reader)
#print(df)
#clippedReadQual(df)
#close(reader)