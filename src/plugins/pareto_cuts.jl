import Plots, Colors


# function pareto_selection_update(V::ConvexApproximation, cut::Cut)
#     # check if cut dominates all other cuts
#     # for (key, coef) in cut.coefficients
#     #     dominates_part = coef >
#     # end

#     for pareto_cut in V.cuts
#         intercept = pareto_cut.intercept
#         cut_coef_dominated = true
#         for (key, coef) in pareto_cut.coefficients
#             cut_coef_dominated = cut_coef_dominated && (cut.coefficients[key] <= coef)
#         end
#         cut_is_dominated = (cut.intercept <= intercept) && cut_coef_dominated 
#         if cut_is_dominated
#             return true
#         end
#     end
#     return false
# end


function pareto_cut_double_pass(V::ConvexApproximation, cut::Cut)
    recalc_min_max_cut_values!(V)
    new_cut_intercept_scaled = scale(
        cut.intercept,
        V.max_cut_values["intercept"],
        V.min_cut_values["intercept"],
    )
    new_cut_coef_scaled = Dict(
        state=> scale(coef, V.max_cut_values["coefs"][state], V.max_cut_values["coefs"][state])
        for state in keys(cut.coefficients)
    )
    idx_cur_used_but_dom_cuts = []
    for (idx, pareto_cut) in enumerate(V.cuts)
        scaled_intercept = scale(
            pareto_cut.intercept,
            V.max_cut_values["intercept"],
            V.min_cut_values["intercept"],
        )
        cut_coef_dominated = true
        for (state, coef_scaled) in new_cut_coef_scaled
            cur_used_cut_coef_scaled = scale(
                pareto_cut.coefficients[state],V.max_cut_values["coefs"][state], V.max_cut_values["coefs"][state]
            )
            # check if even when scaled with epsilon, is new cut coef dominated
            new_cut_coef_dominated = ((1+settings.get("epsilon"))*coef_scaled < cur_used_cut_coef_scaled)
            cut_coef_dominated = cut_coef_dominated && new_cut_coef_dominated
        end
        cut_is_dominated = ((1 + settings.get("epsilon"))*new_cut_intercept_scaled < scaled_intercept) && cut_coef_dominated 
        if cut_is_dominated
            return true
        else
            push!(idx_cur_used_but_dom_cuts, idx)
        end
    end
    cut.pareto_dominant = true
    subproblem = JuMP.owner_model(V.theta)
    for dominated_idx in idx_cur_used_but_dom_cuts
        dominated_cut = V.cuts[dominated_idx]
        JuMP.delete(subproblem, dominated_cut.constraint_ref)
        cut.constraint_ref = nothing
        cut.non_dominated_count = 0
    end
    deleteat!(V.cuts, idx_cur_used_but_dom_cuts)
    # V.cuts = V.cuts[setdiff(1:length(V.cuts), idx_cur_used_but_dom_cuts)]
    # return number of deleted cuts
    return length(idx_cur_used_but_dom_cuts)
end


function cut_is_dominated(V::ConvexApproximation, cut::Cut)::Bool
    new_cut_intercept_scaled = scale(
        cut.intercept,
        V.max_cut_values["intercept"],
        V.min_cut_values["intercept"],
    )
    new_cut_coef_scaled = Dict(
        state=> scale(cut.coefficients[state], V.max_cut_values["coefs"][state], V.max_cut_values["coefs"][state])
        for state in keys(cut.coefficients)
    )
    for pareto_cut in V.cuts
        scaled_intercept = scale(
            pareto_cut.intercept,
            V.max_cut_values["intercept"],
            V.min_cut_values["intercept"],
        )
        cut_coef_dominated = true
        for (state, coef_scaled) in new_cut_coef_scaled
            cur_used_cut_coef_scaled = scale(
                pareto_cut.coefficients[state],V.max_cut_values["coefs"][state], V.max_cut_values["coefs"][state]
            )
            # check if even when scaled with epsilon, is new cut coef dominated
            new_cut_coef_dominated = ((1+settings.get("epsilon"))*coef_scaled < cur_used_cut_coef_scaled)
            cut_coef_dominated = cut_coef_dominated && new_cut_coef_dominated
        end
        cut_is_dominated = ((1 + settings.get("epsilon"))*new_cut_intercept_scaled < scaled_intercept) && cut_coef_dominated 
        if cut_is_dominated
            return true
        end
    end
    cut.pareto_dominant = true
    return false
end

function pareto_frontier!(cuts::Array{Cut})
    # es muss noch ein Argument übergeben werden welches angibt in welcher Dimension die Daten gestackt sind
    # sort data via first entry
    sort!(cuts, by=c->c.intercept, rev=true)

    max_prev_coef = cuts[1].coefficients[:x]
    cuts[1].pareto_dominant = true
    last_pareto = cuts[1]

    for cut in cuts[begin+1:end]
        intercept = cut.intercept
        coef = cut.coefficients[:x]
        if coef > max_prev_coef
            # If next entry has same first entry but higher second entry
            # that entry dominates previous entries
            if intercept == last_pareto.intercept
                last_pareto.pareto_dominant = false
            end
            cut.pareto_dominant = true
            max_prev_coef = coef
            last_pareto = cut
        end
    end
end

function first_cut_dominates(first_cut, second_cut, ValueFunctionApprox)
    cur_min_intercept = ValueFunctionApprox.min_cut_values["intercept"]
    cur_max_intercept = ValueFunctionApprox.max_cut_values["intercept"]
    cur_min_coefs = ValueFunctionApprox.min_cut_values["coefs"]
    cur_max_coefs = ValueFunctionApprox.max_cut_values["coefs"]

    # es wird angenommen das bei gleichen Koeffizienten der aktuellere Cut verwendet wird, da man sonst den gleichen Cut mehrmals im Problem hätte
    for state in keys(first_cut.coefficients)
        scaled_first_cut_coef = scale(first_cut.coefficients[state], cur_max_coefs[state], cur_min_coefs[state])
        scaled_second_cut_coef = scale(second_cut.coefficients[state], cur_max_coefs[state], cur_min_coefs[state])
        # println("epsilon $(settings.get("epsilon"))","First cut",first_cut.coefficients[state], " scaled ", scaled_first_cut_coef, " second cut", second_cut.coefficients[state], " ",  scaled_second_cut_coef)
        # println("dominance relation ",  (((1 + settings.get("epsilon")) * scaled_first_cut_coef) >= scaled_second_cut_coef))
        if ! (((1 + settings.get("epsilon")) * scaled_first_cut_coef) >= scaled_second_cut_coef)
            # first cut does not dominate in all dimensions break for loop
            return false
        end
    end
    scaled_first_intercept = scale(first_cut.intercept, cur_max_intercept, cur_min_intercept) 
    scaled_second_intercept = scale(second_cut.intercept, cur_max_intercept, cur_min_intercept) 

    first_intercept_dominates = (1 + settings.get("epsilon"))*scaled_first_intercept >= scaled_second_intercept
    return first_intercept_dominates

    # frühere Version
    # first_cut_coefs = values(first_cut.coefficients)
    # second_cut_coefs = values(second_cut.coefficients)
    # first_coef_dominates = all(
    #     ((1 + settings.get("epsilon")) .* first_cut_coefs) .>= second_cut_coefs
    # )

    # first_cut_dominates =  (1 + settings.get("epsilon"))*first_cut.intercept >= second_cut.intercept && first_coef_dominates
    # return first_cut_dominates
end

module Bskytree
  import SDDP
  using CxxWrap
  @wrapmodule(
    () -> joinpath("/Users/janik/Documents/Master/KIT/Masterarbeit/src/BSkyTreeJulia/build/lib/","libbskytree_$(SDDP.settings.get("num_dims"))d.dylib")
  )
  function __init__()
    @initcxx()
  end

end # module

"""
    pareto_bskytree(data)

Compute the Pareto front of the given data using the Bskytree algorithm.
Computes only the indices of the elements of the pareto front.
"""
function pareto_bskytree(data::Array{SDDP.Cut})
    # turn data into matrix
    @time matrix = get_cut_point_matrix(data)
    # println("matrix")
    # println(matrix)
    @time pareto_indices = Bskytree.compute_skyline(matrix)
    # println(pareto_indices)
    # println("\\\\\\\\\\\\\\\\\\\\\\")
    return pareto_indices
end


# using Base.Threads
# function get_cut_point_matrix(cuts::Array{SDDP.Cut})::Matrix
#     # Extract intercepts and coefficients directly
#     n_cuts = length(cuts)
#     n_coeffs = length(first(cuts).coefficients)
#     matrix = zeros(n_cuts, n_coeffs + 1)

#     @threads for i in 1:n_cuts
#         matrix[i, 1:n_coeffs] .= collect(values(cuts[i].coefficients))
#         matrix[i, end] = cuts[i].intercept
#     end

#     return matrix

#     return matrix
# end

function get_cut_point_matrix(cuts::Array{SDDP.Cut})::Matrix
    cut_coefs = values.(getfield.(cuts, :coefficients))
    matrix = zeros(length(cuts), length(cut_coefs[1])+1)

    for (row, coefs) in enumerate(cut_coefs)
        for (col, coef) in enumerate(coefs)
            matrix[row, col] = coef
            matrix[row, end] = cuts[row].intercept
        end
    end
    return matrix
end


"""
    bnl!(data)

Compute the Pareto front of the given data using the Block Nested Loop (BNL) algorithm.

# Arguments
- `data`: An array of cuts representing the data.

# Returns
An array of cuts representing the Pareto front.

"""
function bnl!(ValueFunctionApprox::ConvexApproximation)
    cut_array = filter(cut -> cut.constraint_ref !== nothing, ValueFunctionApprox.cuts)
    if isempty(cut_array)
        return cut_array
    elseif length(cut_array) == 1
        cut_array[1].pareto_dominant = true
        return cut_array
    end
    window = [cut_array[1]]
    window[1].pareto_dominant = true

    println("=====================")
    for idx in 2:length(cut_array)
        cut = cut_array[idx]

        dominated = false
        i = 1
        while i <= length(window)
            window_cut = window[i]
            # check if a window cut dominates the new cut, if so discard cut
            if first_cut_dominates(window_cut, cut, ValueFunctionApprox)
                dominated = true
                cut.pareto_dominant = false
                break
            # check if the new cut dominates a the window cut, if so delete window cut
            # otherwise cut can be added to window (aka. block)
            elseif first_cut_dominates(cut, window_cut, ValueFunctionApprox)
                println("dominated cut in window")
                window_cut.pareto_dominant = false
                # window[i].pareto_dominant = false
                deleteat!(window, i)
            else
                i+=1
            end
        end
        if !dominated
            cut.pareto_dominant = true
            push!(window, cut)
        end
    end

    println("Deleted Cuts in bnl ",length(ValueFunctionApprox.cuts) - length(window), " individual length $(length(ValueFunctionApprox.cuts)), $(length(window))")
    for cut in window
        if !cut.pareto_dominant
            println("not dominant")
        end
        println(cut.constraint_ref)
    end
    return window
end

function recalc_min_max_cut_values!(ValueFunctionApprox::ConvexApproximation)
    # Recalc min max value after pareto cuts were recalculated
    if isempty(ValueFunctionApprox.cuts)
        ValueFunctionApprox.min_cut_values["intercept"] = Inf
        ValueFunctionApprox.max_cut_values["intercept"] = -Inf
        for state in keys(ValueFunctionApprox.states)
            ValueFunctionApprox.min_cut_values["coefs"][state] = Inf
            ValueFunctionApprox.max_cut_values["coefs"][state] = -Inf
        end

        return
    end
    first_cut = ValueFunctionApprox.cuts[1]
    for state in keys(first_cut.coefficients)
        ValueFunctionApprox.min_cut_values["coefs"][state] = first_cut.coefficients[state]
        ValueFunctionApprox.max_cut_values["coefs"][state] = first_cut.coefficients[state]
    end
    ValueFunctionApprox.min_cut_values["intercept"] = first_cut.intercept
    ValueFunctionApprox.max_cut_values["intercept"] = first_cut.intercept

    for cut in ValueFunctionApprox.cuts[2:end]
        # Optionally mark as dominant (if not already)
        cut.pareto_dominant = true
    
        for state in keys(cut.coefficients)
            ValueFunctionApprox.min_cut_values["coefs"][state] = min(ValueFunctionApprox.min_cut_values["coefs"][state], cut.coefficients[state])
            ValueFunctionApprox.max_cut_values["coefs"][state] = max(ValueFunctionApprox.max_cut_values["coefs"][state], cut.coefficients[state])
        end
    
        ValueFunctionApprox.min_cut_values["intercept"] = min(ValueFunctionApprox.min_cut_values["intercept"], cut.intercept)
        ValueFunctionApprox.max_cut_values["intercept"] = max(ValueFunctionApprox.max_cut_values["intercept"], cut.intercept)
    end 
end



scale(x, max_val, min_val) = (x-min_val)/(max_val - min_val)
"""
Heursitically check if a point is dominated. That is checked by scaling the intercepts and coefficients of the array to 
the interval [0, 1] and check if a point falls below the 75% lines for every dimension. The point is then considered to be
dominated, as better points will typically lay in 
"""
function heuristic_point_is_dominated(cut::Cut, V::ConvexApproximation ; threshold::Float64)::Bool
    # check intercept
    scaled_intercept = scale(cut.intercept, V.max_cut_values["intercept"], V.min_cut_values["intercept"])
    if scaled_intercept >= threshold
        return false
    end
    # check coefs
    for state in keys(V.states)
        scaled_coef = scale(
            cut.coefficients[state],
            V.max_cut_values["coefs"][state],
            V.min_cut_values["coefs"][state],
        )
        if scaled_coef >= threshold
            return false
        end
    end
    return true
end



function viz_pareto_cuts(cuts::Array{Cut})::Plots.Plot
    color_palette = Colors.distinguishable_colors(length(cuts))
    cut_dominates = getfield.(cuts, :pareto_dominant)
    marker_mask = ifelse.(cut_dominates, :star6, :circle)
    pareto_sicht = Plots.scatter(
        map(c->c.intercept, cuts),
        map(c->c.coefficients[:x], cuts),
        color=color_palette,
        title="Slope intercept",
        ylabel="Intercept",
        xlabel="Slope",
        legend=false,
        markersize=10,
        markershape=marker_mask,
    )
    linear_cuts = Plots.plot(title="Cuts", legend=false)
    for (i,cut) in enumerate(cuts)
        linear_cut(x) = cut.intercept + cut.coefficients[:x]*x
        Plots.plot!(linear_cuts, 0:0.01:.25, linear_cut, color=color_palette[i], linewidth=3)
    end 
    Plots.plot!([0], seriestype=:vline, color=:red)

    p = Plots.plot(linear_cuts, pareto_sicht, layout=(1,2), size=(1000,500), linewidth=2)
    return p 
end

