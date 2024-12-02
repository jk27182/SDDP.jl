function _write_cut_buffer(file_name::String)
    open(file_name, "w") do file
        return write(file, CUT_BUFFER)
    end
end

function _cut_selection_update(
    V::ConvexApproximation,
    cut::Cut,
    state::Dict{Symbol,Float64};
)
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
        # println("Cuts werden geloescht $(length(V.cuts_to_be_deleted))")
        # println("delete $(length(V.cuts_to_be_deleted)) cuts")
        # SDDP.cut_deletion_standard_per_stage[stage] += length(V.cuts_to_be_deleted)
        for cut in V.cuts_to_be_deleted
            JuMP.delete(model, cut.constraint_ref)
            cut.constraint_ref = nothing
            cut.non_dominated_count = 0
        end
    end
    empty!(V.cuts_to_be_deleted)
    return
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
        cut1.intercept == cut2.intercept &&
        cut1.coefficients == cut2.coefficients &&
        cut1.obj_y == cut2.obj_y &&
        cut1.belief_y == cut2.belief_y &&
        cut1.non_dominated_count == cut2.non_dominated_count &&
        cut1.constraint_ref == cut2.constraint_ref
    )
end


"""
    prune_cuts!(model::PolicyGraph)

    This function iterates over all nodes in the model and removes the cuts that are not Pareto dominant.
    A cut is considered Pareto dominant if there is no other cut that is better in both the coefficients and intercept.
"""
# function prune_cuts!(::SDDP.SINGLE_CUT, model::PolicyGraph{T}) where T
#     Logging.@debug("ich bin in prune cuts!!!!!!!!!!!")
#     Logging.@debug("-----------------")
#     for (stage, node) in model.nodes
#         Logging.@debug("-----------------")
#         prune_cuts_inner!(node.bellman_function.global_theta)
#         # cuts = node.bellman_function.global_theta.cuts
#         # pareto_dominant_cuts = SDDP.bnl!(cuts)
#         # # owner_model(model[1].bellman_function.global_theta.cuts[6].constraint_ref)
        
#         # Logging.@debug("In stage $stage  #Cuts that need to be deleted: $(length(filter(cut -> !cut.pareto_dominant, cuts)))")
        
#         # for dominated_cut in filter(cut -> !cut.pareto_dominant, cuts)
#         #     Logging.@debug("Cut that needs to be deleted: $dominated_cut")
#         #     if dominated_cut.constraint_ref !== nothing
#         #         JuMP.delete(node.subproblem, dominated_cut.constraint_ref)
#         #         dominated_cut.constraint_ref = nothing
#         #         dominated_cut.non_dominated_count = 0
#         #     end
#         # end
#         # node.bellman_function.global_theta.cuts = pareto_dominant_cuts
#     end
# end



function prune_cuts!(model::PolicyGraph{T}) where T
    # prune_cuts_inner!(approx) = prune_cuts_inner_h!(approx)
    pruning_type = settings.get("pruning_type") 
    if pruning_type == "bnl"
        prune_cuts_inner! = prune_cuts_inner_bnl!
    elseif pruning_type == "Bskytree"
        prune_cuts_inner! = prune_cuts_inner_bskytree!
    else
        prune_cuts_inner! = prune_cuts_inner_heuristic!
    end

    n_stages = length(model.nodes)
    for (stage, node) in model.nodes
        if stage == n_stages || stage  == string(2)
            continue
        end
        node = model.nodes[stage]
        cut_type = node.bellman_function.cut_type
        if cut_type == SDDP.MULTI_CUT
            # loops through all scenarios for the respective approximation
            for local_approx in node.bellman_function.local_thetas
                prune_cuts_inner!(local_approx)
            end
        else
            # else assume Single Cuts
            global_approx = node.bellman_function.global_theta
            prune_cuts_inner!(global_approx)
            # prune_cuts_inner_h!(global_approx)
        end
    end
    # Logging.@debug("-----------------")
end


"""
remove element at index idx from arr in O(1) by ignoring ordering in they array.
"""
function remove_at_idx_unordered!(arr::Array, idx)::Array
    if !(1<= idx <= length(arr))
        throw(BoundsError("Index $(idx) out of bounds [1, $(length(arr))]"))
    end
    if length(arr) > 1
        arr[idx], arr[end] = arr[end], arr[idx]
    end
    pop!(arr)
end

function prune_cuts_inner_bskytree!(ValueFunctionApprox::ConvexApproximation)
    println("RUNNING Bskytree")
    pareto_indices = SDDP.pareto_bskytree(
        ValueFunctionApprox.cuts
        # filter(cut -> cut.constraint_ref !== nothing, ValueFunctionApprox.cuts)
    )
    deletion_idxs = setdiff(1:length(ValueFunctionApprox.cuts), pareto_indices)
    subproblem = JuMP.owner_model(ValueFunctionApprox.theta)
    @time begin
        for del_idx in deletion_idxs
            JuMP.delete(subproblem, ValueFunctionApprox.cuts[del_idx].constraint_ref)
        end
    end
    @time deleteat!(ValueFunctionApprox.cuts, deletion_idxs)

    # for p_idx in pareto_indices
    #     ValueFunctionApprox.cuts[p_idx].pareto_dominant = true
    # end
    # subproblem = JuMP.owner_model(ValueFunctionApprox.theta)
    # for (idx, cut) in enumerate(ValueFunctionApprox.cuts)
    #     # Logging.@debug("Cut that needs to be deleted: $dominated_cut")
    #     if cut.pareto_dominant == false
    #         println("delete cut")
    #         println(cut.constraint_ref)
    #         JuMP.delete(subproblem, cut.constraint_ref)
    #         # remove_at_idx_unordered!(ValueFunctionApprox.cuts, idx)
    #         cut.constraint_ref = nothing
    #         cut.non_dominated_count = 0
    #     end
    # end
    # ValueFunctionApprox.cuts = ValueFunctionApprox.cuts[pareto_indices]
end

function prune_cuts_inner_bnl!(ValueFunctionApprox::ConvexApproximation)
    # Cuts with Constraint ref of nothing are not in the model anymore
    # Kommt nicht vor wenn die Selection Operationen passend ausgewählt werden
    println("RUNNING BNL")
    pareto_dominant_cuts = SDDP.bnl!(
        ValueFunctionApprox,
    )

    # if length(cuts_to_delete) > 0
    #     # println("$(length(cuts_to_delete)) cut(s) need to be deleted")
    # end 
    println(repeat("=", 10))
    # println("# cuts", length(ValueFunctionApprox.cuts))
    # println("# pareto dominant cuts", length(pareto_dominant_cuts))
    # for cut in ValueFunctionApprox.cuts
    #     println(cut.constraint_ref)
    # end
    # println(repeat("=", 10), "Pareto Dominant Cuts")
    # for cut in pareto_dominant_cuts
    #     println(cut.constraint_ref)
    # end
    # println(repeat("=", 10))
    println("start deleting stuff")
    for (idx, cut) in enumerate(ValueFunctionApprox.cuts)
        # Logging.@debug("Cut that needs to be deleted: $dominated_cut")
        if cut.pareto_dominant == false
            println("delete cut")
            println(cut.constraint_ref)
            subproblem = JuMP.owner_model(ValueFunctionApprox.theta)
            JuMP.delete(subproblem, cut.constraint_ref)
            # remove_at_idx_unordered!(ValueFunctionApprox.cuts, idx)
            cut.constraint_ref = nothing
            cut.non_dominated_count = 0
        end
    end
    println("end deleting stuff")
    println(repeat("=", 10))
    ValueFunctionApprox.cuts = pareto_dominant_cuts
end

function prune_cuts_inner_heuristic!(ValueFunctionApprox::ConvexApproximation)
    cur_min_intercept = ValueFunctionApprox.min_cut_values["intercept"]
    cur_max_intercept = ValueFunctionApprox.max_cut_values["intercept"]
    cur_min_coefs = ValueFunctionApprox.min_cut_values["coefs"]
    cur_max_coefs = ValueFunctionApprox.max_cut_values["coefs"]

    # cur_min_coefs = Dict()
    # cur_max_coefs = Dict()
    # for state in keys(ValueFunctionApprox.states)
    #     cur_min_coefs[state] = Inf
    #     cur_max_coefs[state] = -Inf
    # end
    # cur_max_intercept = -Inf
    # cur_min_intercept = Inf
    for cut in ValueFunctionApprox.cuts
        # sollte min und max vals nach jeder iteration angepasst werden?
        # mache es bewusst erstmal nicht
        # muss nochmal drüber nachdenken, ob sonst zu viele Punkte gelöscht werden
        # ist das überhaupt möglich?
        if cut.constraint_ref !== nothing
            is_dominated = SDDP.heuristic_point_is_dominated(
                cut,
                ValueFunctionApprox,
                threshold=settings.get("threshold"),
            )
            if is_dominated
                # println("delete cut")
                # print("point is dominated")
                # is_dominated = SDDP.heuristic_point_is_dominated(
                #     cut,
                #     ValueFunctionApprox,
                #     threshold=settings.get("threshold"),
                # )
                cut.pareto_dominant = false
                # if cut is deleted, minimum and maximum need to be updated
                subproblem = JuMP.owner_model(ValueFunctionApprox.theta)
                JuMP.delete(subproblem, cut.constraint_ref)
                cut.constraint_ref = nothing
                cut.non_dominated_count = 0
            else
                # ist falsch denke ich, so wird das nicht funktionieren

            
                # if cut is not dominated, hence not removed, check if cut is new min or max for subsequent pruning steps
                cur_min_intercept = min(cut.intercept, cur_min_intercept)
                cur_max_intercept = max(cut.intercept, cur_max_intercept)
                for state in keys(ValueFunctionApprox.states)
                    cur_min_coefs[state] = min(cut.coefficients[state], cur_min_coefs[state])
                    cur_max_coefs[state] = max(cut.coefficients[state], cur_max_coefs[state])
                end
            end
        end
    end
    # cut_matrix = SDDP.get_cut_point_matrix(filter(cut-> cut.constraint_ref !== nothing, ValueFunctionApprox.cuts))
    # ValueFunctionApprox.min_cut_values = minimum(cut_matrix, dim=1)
    # ValueFunctionApprox.max_cut_values = maximum(cut_matrix, dim=1)
    ValueFunctionApprox.min_cut_values["intercept"] = cur_min_intercept
    ValueFunctionApprox.max_cut_values["intercept"] = cur_max_intercept
    ValueFunctionApprox.min_cut_values["coefs"] = cur_min_coefs
    ValueFunctionApprox.max_cut_values["coefs"] = cur_max_coefs
end
