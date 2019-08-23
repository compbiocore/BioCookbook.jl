module BioCookbook

import FASTX
import BioSequences

export concatenate, simAlignments

include("concatenate.jl")
include("simulation.jl")

end