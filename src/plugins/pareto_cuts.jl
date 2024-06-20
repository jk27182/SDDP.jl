# using Plots

function cut_is_dominated(V::ConvexApproximation, cut::Cut)
    # pareto_front = V.cuts
    for pareto_cut in V.cuts
        intercept = pareto_cut.intercept
        # use value of dictionary for comparison in all
        # check coef domination
        # assumes that atleast one coefficient exist, otherwise 
        cut_coef_dominated = true
        for (key, coef) in pareto_cut.coefficients
            cut_coef_dominated = cut_coef_dominated && (cut.coefficients[key] <= coef)
        end
        # coef_dominated = all(
        #     pair -> pair[2] <= coef,
        #     pareto_cut.coefficients
        # )
        cut_is_dominated = (cut.intercept <= intercept) && cut_coef_dominated 
        # coef = pareto_cut.coefficients[:x_storage]
        # cut_is_dominated = (cut.intercept <= intercept) && (cut.coefficients[:x_storage]  <= coef)
        if cut_is_dominated
            return true
        end
    end
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

function viz_pareto_cuts(cuts::Array{Cut})::Plots.Plot
    color_palette = distinguishable_colors(length(cuts))
    cut_dominates = getfield.(cuts, :pareto_dominant)
    marker_mask = ifelse.(cut_dominates, :star6, :circle)
    pareto_sicht = scatter(
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
    linear_cuts = plot(title="Cuts", legend=false)
    for (i,cut) in enumerate(cuts)
        linear_cut(x) = cut.intercept + cut.coefficients[:x]*x
        plot!(linear_cuts, 0:0.01:.25, linear_cut, color=color_palette[i], linewidth=3)
    end 
    plot!([0], seriestype=:vline, color=:red)

    p = plot(linear_cuts, pareto_sicht, layout=(1,2), size=(1000,500), linewidth=2)
    return p 
end

