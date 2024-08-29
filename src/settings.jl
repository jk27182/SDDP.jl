mutable struct Settings
    use_pareto_cut_logic::Bool
    log_level::Int
    use_cut_selection::Bool
    use_pruning::Bool
    prune_interval::Int

    # use a closure to store the settings and prevent public acces
    function Settings(
        ;
        use_pareto_cut_logic::Bool,
        log_level::Int,
        use_cut_selection::Bool,
        use_pruning::Bool,
        prune_interval::Int,
    )
        settings = new(
            use_pareto_cut_logic,
            log_level,
            use_cut_selection,
            use_pruning,
            prune_interval
        )
        function set(;use_pareto_cut_logic::Bool, log_level::Int, use_cut_selection::Bool, use_pruning::Bool, prune_interval::Int)
            settings.use_pareto_cut_logic = use_pareto_cut_logic
            settings.log_level = log_level
            settings.use_cut_selection = use_cut_selection
            settings.use_pruning = use_pruning
            settings.prune_interval = prune_interval
            return
        end
        set(field::Symbol, value) = setfield!(settings, field, value)
        set(field::String, value) = setfield!(settings, Symbol(field), value)

        get() = settings
        get(x::Symbol) = getfield(settings, x)
        get(x::String) = getfield(settings, Symbol(x))

        return (;set, get)
    end
end
