function _write_cut_buffer(file_name::String)
    open(file_name, "w") do file
        write(file, CUT_BUFFER)
    end 
end

function _add_cut(
    V::ConvexApproximation,
    θᵏ::Float64, #intercept
    πᵏ::Dict{Symbol,Float64},
    xᵏ::Dict{Symbol,Float64}, #state
    obj_y::Union{Nothing,NTuple{N,Float64}},
    belief_y::Union{Nothing,Dict{T,Float64}};
    cut_selection::Bool = true,
    cut_buffering::Bool,
) where {N,T}
    for (key, x) in xᵏ
        θᵏ -= πᵏ[key] * x
    end
    println("Das laeuft in", @__FILE__())
    _dynamic_range_warning(θᵏ, πᵏ)
    cut = Cut(θᵏ, πᵏ, obj_y, belief_y, 1, nothing)
    push!(CUT_BUFFER, cut)
    if !cut_buffering 
        _add_cut_constraint_to_model(V, cut)
        if cut_selection
            _cut_selection_update(V, cut, xᵏ)
        end
    end
    # _add_cut_constraint_to_model(V, cut)

    push!(CUT_BUFFER, cut)
    if cut_selection
        _cut_selection_update(V, cut, xᵏ)
    end


    return
end
