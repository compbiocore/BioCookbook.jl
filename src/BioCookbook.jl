module BioCookbook

import BioSequences
import BioSequences
import BioAlignments
import Statistics

module msa
import FASTX
include("concatenate.jl")
include("simulation.jl")
export concatenate, simAlignments
end

using .msa

module mapping
import XAM
include("extractReads.jl")
export extractReads

module viz
import XAM
import DataFrames
include("clipanalysis.jl")
end

end

using .mapping

end