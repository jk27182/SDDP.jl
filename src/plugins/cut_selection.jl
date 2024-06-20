function _write_cut_buffer(file_name::String)
    open(file_name, "w") do file
        write(file, CUT_BUFFER)
    end 
end


function _cut_selection_update(
    V::ConvexApproximation,
    cut::Cut,
    state::Dict{Symbol,Float64};
    iteration::Int,
    stage::Int,
)
    SDDP.count_update += 1
    # println("hier laeuft die cut selection")
    model = JuMP.owner_model(V.theta)
    is_minimization = JuMP.objective_sense(model) == MOI.MIN_SENSE
    # Ein sampled State entspricht dem x mit dem dann ein Cut generiert wird
    # println("der belief ist $(cut.belief_y)")
    # belief state war in Newsvendor nothin, macht auch Sinn da der Belief State etwas mit partially observable States zu tun hat
    # In Newsvendor ist der State fully observable
    sampled_state = SampledState(state, cut.obj_y, cut.belief_y, cut, NaN)
    sampled_state.best_objective = _eval_height(cut, sampled_state)
    # println("Die bisherige Approximation ist $(V)")
    # println("Der SampledState ist $(sampled_state)")
    # Loop through previously sampled states and compare the height of the most
    # recent cut against the current best. If this new cut is an improvement,
    # store this one instead.
    for old_state in V.sampled_states
        # Only compute cut selection at same points in concave space.
        if old_state.obj_y != cut.obj_y || old_state.belief_y != cut.belief_y
            continue
        end
        height = _eval_height(cut, old_state)
        if _dominates(height, old_state.best_objective, is_minimization)
            old_state.dominating_cut.non_dominated_count -= 1
            cut.non_dominated_count += 1
            old_state.dominating_cut = cut
            old_state.best_objective = height
        end
    end
    push!(V.sampled_states, sampled_state)
    # Now loop through previously discovered cuts and compare their height at
    # `sampled_state`. If a cut is an improvement, add it to a queue to be
    # added.
    for old_cut in V.cuts
        if old_cut.constraint_ref !== nothing
            # We only care about cuts not currently in the model.
            continue
        elseif old_cut.obj_y != sampled_state.obj_y
            # Only compute cut selection at same points in objective space.
            continue
        elseif old_cut.belief_y != sampled_state.belief_y
            # Only compute cut selection at same points in belief space.
            continue
        end
        height = _eval_height(old_cut, sampled_state)
        if _dominates(height, sampled_state.best_objective, is_minimization)
            sampled_state.dominating_cut.non_dominated_count -= 1
            old_cut.non_dominated_count += 1
            sampled_state.dominating_cut = old_cut
            sampled_state.best_objective = height
            _add_cut_constraint_to_model(V, old_cut)
        end
    end
    # push!(CUT_DICT[stage], cut)
    push!(V.cuts, cut)
    # Delete cuts that need to be deleted.
    for cut in V.cuts
        if cut.non_dominated_count < 1
            if cut.constraint_ref !== nothing
                push!(V.cuts_to_be_deleted, cut)
            end
        end
    end
    if length(V.cuts_to_be_deleted) >= V.deletion_minimum
        for cut in V.cuts_to_be_deleted
            JuMP.delete(model, cut.constraint_ref)
            cut.constraint_ref = nothing
            cut.non_dominated_count = 0
        end
    end
    empty!(V.cuts_to_be_deleted)
    return
end

function _cut_selection_update_pareto_efficient(
    V::ConvexApproximation,
    cut::Cut,
    state::Dict{Symbol,Float64},
    iteration::Int;
    stage::Int,
)
    println("neue cut selection mit multiple dispatch")
    SDDP.count_update += 1
    # println("hier laeuft die cut selection")
    if !cut_is_dominated(V.cuts, cut)
        cut.pareto_dominant = true
        push!(V.cuts, cut)
    end
    a = 1
end

function _cut_selection_update_pareto(
    V::ConvexApproximation,
    cut::Cut,
    state::Dict{Symbol,Float64},
    iteration::Int;
    stage::Int,
)
    println("neue cut selection mit multiple dispatch")
    SDDP.count_update += 1
    # println("hier laeuft die cut selection")
    push!(V.cuts, cut)
    cuts_copy = deepcopy(V.cuts)
    pareto_frontier!(cuts_copy)
    filter!(c->c.pareto_dominant, cuts_copy) 
    a = 1
end

# In V sind alle bisherigen Cuts abgepseichert und dann wird das ganze erweitert
# mit dem neuen Cut 'cut' und dem momentanen State 'state'
# function _cut_selection_update(
#     V::ConvexApproximation,
#     cut::Cut,
#     state::Dict{Symbol,Float64},
#     iteration::Int;
#     stage::Int,
# )
#     SDDP.count_update += 1
#     model = JuMP.owner_model(V.theta)
#     is_minimization = JuMP.objective_sense(model) == MOI.MIN_SENSE
#     # Ein sampled State entspricht dem x mit dem dann ein Cut generiert wird
#     # println("der belief ist $(cut.belief_y)")
#     # belief state war in Newsvendor nothin, macht auch Sinn da der Belief State etwas mit partially observable States zu tun hat
#     # In Newsvendor ist der State fully observable
#     sampled_state = SampledState(state, cut.obj_y, cut.belief_y, cut, NaN)
#     sampled_state.best_objective = _eval_height(cut, sampled_state)
#     # test_cut = SDDP.Cut(0.0, Dict(:x => 1.0), nothing, nothing, 1, nothing)

#     push!(V.cuts, cut)
#     for cut in V.cuts
#         if !cut_is_in(cut, pareto_cuts)
#             println("Der Cut $(cut.constraint_ref) wird geloescht")
#             JuMP.delete(model, cut.constraint_ref)
#             # if cut.constraint_ref !== nothing
#             # end
#         end
#     end
# end

function cut_is_in(cut::Cut, pareto_cuts::Array{Cut})
    for pareto_cut in pareto_cuts
        if is_equal_cut(cut, pareto_cut)
            return true
        end
    end
    return false
end


function is_equal_cut(cut1::Cut, cut2::Cut)
    return (
        cut1.intercept == cut2.intercept 
        && cut1.coefficients == cut2.coefficients
        && cut1.obj_y == cut2.obj_y
        && cut1.belief_y == cut2.belief_y 
        && cut1.non_dominated_count == cut2.non_dominated_count 
        && cut1.constraint_ref == cut2.constraint_ref 
    )
end
