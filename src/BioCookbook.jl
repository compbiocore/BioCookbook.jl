module BioCookbook

import FASTX
import XAM
import BioSequences
import BioAlignments
import Statistics
import DataFrames

export concatenate, simAlignments,
	extractReads

include("concatenate.jl")
include("simulation.jl")
include("extractReads.jl")
include("clipanalysis.jl")

end