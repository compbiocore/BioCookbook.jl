function addRecords!(dest::FASTX.FASTA.Record, src::FASTX.FASTA.Record)
  append!(dest.data, src.data[src.sequence])
  dest.filled = UnitRange{Int}(dest.filled.start, dest.filled.stop + src.filled.stop - src.sequence.start + 1)
  dest.sequence = UnitRange{Int}(dest.sequence.start, dest.sequence.stop + src.sequence.stop - src.sequence.start + 1)
end

function concatenate(infiles::Vector{String}, outfile::String)
  concat = Dict{String, FASTX.FASTA.Record}()
  allLabels = Set(String[])
  alignLen = 0
  alphabet=""
  firstrecord=true
  for (i, file) in enumerate(infiles)
    theseLabels = String[]
    if i == 1 # first file will have only new ids
      reader = open(FASTX.FASTA.Reader, file)
      for record in reader
        if firstrecord==true
          alignLen += length(record.sequence)
          firstrecord=false
          alphabet = typeof(FASTX.FASTA.sequence(record))
        end
        id = FASTX.FASTA.identifier(record)
        push!(allLabels, id)
        concat[id] = record
      end
      close(reader)
    else # later records need to check for id existance
      record = FASTX.FASTA.Record()
      open(FASTX.FASTA.Reader, file) do reader
        while !eof(reader)
          read!(reader, record)
          id = FASTX.FASTA.identifier(record)
          push!(theseLabels, id)
          if id in allLabels
            addRecords!(concat[id], record)
          else # need to create and back-propogate gaps
            push!(allLabels, id)
            newRecord = FASTX.FASTA.Record(id,repeat("-",alignLen))
            pop!(newRecord.data)
            concat[id] = newRecord
            addRecords!(concat[id], record)
          end
        end
        alignLen += length(record.sequence)
        # assume last sequence will be representative of alphabet
        @assert alphabet == typeof(FASTX.FASTA.sequence(record))
      end
      missingLabels = setdiff(allLabels, theseLabels)
      for label in missingLabels
        gapRecord = copy(record)
        gapRecord.data[gapRecord.sequence] = repeat(UInt8[0x2d], length(gapRecord.sequence))
        addRecords!(concat[label], gapRecord)
      end
    end
  end
  open(FASTX.FASTA.Writer, outfile) do w
    for (key, value) in concat
      write(w, value)
    end
  end
  return concat
end

#using BenchmarkTools
#tmpdir = tempdir()
#infiles = simAlignments(30,100,tmpdir)
#@benchmark concatenate($simAlignments(30,100,tmpdir),$mktemp(tmpdir)[1])
#rm(tmpdir, recursive=true)