#  Copyright (c) 2017-23, Oscar Dowson and SDDP.jl contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

module SDDP

import Reexport
Reexport.@reexport using JuMP

import Distributed
import HTTP
import JSON
import MutableArithmetics
import Printf
import Random
import SHA
import Statistics
import TimerOutputs
import Plots
import Logging

# Work-around for https://github.com/JuliaPlots/RecipesBase.jl/pull/55
# Change this back to `import RecipesBase` once the fix is tagged.
using RecipesBase

export @stageobjective


include("settings.jl")
const settings = Settings(
    use_pareto_cut_logic = true,
    log_level = 1,
    use_cut_selection = true,
    use_pruning = false,
    prune_interval = 10, 
)
# SETTINGS = JSON.Parser.parsefile("/Users/janik/Documents/Master/KIT/Masterarbeit/src/SDDP.jl/src/plugins/settings.json")

logger = Logging.ConsoleLogger(stdout, settings.get("log_level"))
Logging.global_logger(logger)


# Modelling interface.
include("user_interface.jl")
include("modeling_aids.jl")

# Default definitions for SDDP related modular utilities.
include("plugins/headers.jl")

# Tools for overloading JuMP functions
include("binary_expansion.jl")
include("JuMP.jl")

# Printing utilities.
include("cyclic.jl")
include("print.jl")

# The core SDDP code.
include("algorithm.jl")

# Specific plugins.
include("plugins/risk_measures.jl")
include("plugins/sampling_schemes.jl")
include("plugins/bellman_functions.jl")
include("plugins/pareto_cuts.jl")
include("plugins/cut_selection.jl")
include("plugins/stopping_rules.jl")
include("plugins/local_improvement_search.jl")
include("plugins/duality_handlers.jl")
include("plugins/parallel_schemes.jl")
include("plugins/backward_sampling_schemes.jl")
include("plugins/forward_passes.jl")

CUT_LISTE::Dict{Int, Any} = Dict()
CUT_DICT = Dict{Any, Array{Cut}}(node => Cut[] for node in 1:2)
count_update::Int = 0
iter_count::Int = 0

# Visualization related code.
include("visualization/publication_plot.jl")
include("visualization/spaghetti_plot.jl")
include("visualization/dashboard.jl")
include("visualization/value_functions.jl")

# Other solvers.
include("deterministic_equivalent.jl")
include("biobjective.jl")
include("alternative_forward.jl")

include("Experimental.jl")
include("MSPFormat.jl")

end
