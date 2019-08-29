using FASTX
using CodecZlib
# Zlib could possibly be made more efficient with transcodingstreams
# https://github.com/bicycle1885/TranscodingStreams.jl
# NOTE: if you assume that each FASTQ files are sorted in same order then
# you do not have to hash
# NOTE: Implementation is for single set of files
# Assumes file names will be:
# $input_1.fq.gz
# $input_2.fq.gz
# $input_1.tagged_filter.fastq.gz
# $input_2.tagged_filter.fastq.gz
# the search function is used to build the union and map_search is the intersection
function match_screens(input, out; merge_type = "intersection", raw_fq=".fq.gz", filt_fq = ".tagged_filter.fastq.gz", out_fq=".fs.fq.gz")
    r1_map = make_map(string(input,"_1", filt_fq))
    r2_map = make_map(string(input,"_2", filt_fq))
    writer_r1 = FASTQ.Writer(GzipCompressorStream(open(string(out,"_1",out_fq),"w")))
    writer_r2 = FASTQ.Writer(GzipCompressorStream(open(string(out,"_2", out_fq),"w")))
    map_search!(r1_map,r2_map,writer_r1,writer_r2)
    if merge_type == "union"
        search!(string(input,"_2", raw_fq),r1_map,writer_r2)
        search!(string(input,"_1". raw_fq),r2_map,writer_r1)
    end
    close(writer_r1)
    close(writer_r2)
end
# search thru r1_map for records in r2_map. If matched, write to file and
# delete record from both r1_map and r2_map
function map_search!(r1_map, r2_map, w1, w2)
    #i = 0 # TESTING THAT ALL FUNCTIONS ARE WORKING
    for record in values(r2_map)
        header = FASTQ.identifier(record)
        #p = parse(header)
        paired_rec = pop!(r1_map,header,0)
        if paired_rec != 0
            write(w1, record)
            write(w2, paired_rec)
            delete!(r2_map,header)
        end
    #    i = i + 1
    #    if i == 10000
    #        break
    #    end
    end
end
# search thru reader for records with identifiers in rmap
# rmap will have records deleted as we go based on if they are found or not
# For primer on mutable functions in Julia: https://docs.julialang.org/en/latest/manual/faq/?highlight=Mutable#functions
function search!(file, rmap, w1, w2)
    #i = 0
    reader = FASTQ.Reader(GzipDecompressorStream(open(file,"r")))
    for record in reader
        header = FASTQ.identifier(record)
        #p = parse(header)
        paired_rec = pop!(rmap,header,0)
        if paired_rec != 0
            write(w1, record)
            write(w2, paired_rec)
        end
    #    i = i + 1
    #    if i == 10000
    #        break
    #    end
    end
    close(reader)
end
# make read map
# reader is FASTQ Reader object
function make_map(file)
    reader = FASTQ.Reader(GzipDecompressorStream(open(file,"r")))
    rmap = Dict()
    for record in reader
        header = FASTQ.identifier(record)
        #p = parse(header)
        rmap[header] = record
    end
    return(rmap)
end
# parse header in some way that is mappable to the other read
# appears to not be needed for now
#function parse(header)
#end
@time match_screens("/Users/aguang/CORE/scratch/fastq_screen/D1121A","/Users/aguang/CORE/scratch/fastq_screen/test.fq.gz")