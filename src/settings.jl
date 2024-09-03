mutable struct Settings
    use_pareto_cut_logic::Bool
    log_level::Int
    use_cut_selection::Bool
    use_pruning::Bool
    prune_interval::Int
    debug_mode::Bool

    # use a closure to store the settings and prevent public acces
    function Settings(;
        use_pareto_cut_logic::Bool,
        log_level::Int,
        use_cut_selection::Bool,
        use_pruning::Bool,
        prune_interval::Int,
        debug_mode::Bool,
    )
        settings = new(
            use_pareto_cut_logic,
            log_level,
            use_cut_selection,
            use_pruning,
            prune_interval,
            debug_mode,
        )
        function set!(;
            use_pareto_cut_logic::Bool,
            log_level::Int,
            use_cut_selection::Bool,
            use_pruning::Bool,
            prune_interval::Int,
            debug_mode::Bool,
        )
            settings.use_pareto_cut_logic = use_pareto_cut_logic
            settings.log_level = log_level
            settings.use_cut_selection = use_cut_selection
            settings.use_pruning = use_pruning
            settings.prune_interval = prune_interval
            settings.debug_mode = debug_mode
            return
        end
        set!(field::Symbol, value) = setfield!(settings, field, value)
        set!(field::String, value) = setfield!(settings, Symbol(field), value)

        get() = settings
        get(x::Symbol) = getfield(settings, x)
        get(x::String) = getfield(settings, Symbol(x))

        function get_setting_id()
            println("get_setting_id")
            @show settings
            pareto_str = get(:use_pareto_cut_logic) ? "Pareto Cuts" : ""
            cut_selection_str =
                get(:use_cut_selection) ? "Standard Cut Selection" : ""
            pruning_str =
                get(:use_pruning) ?
                "Pruning every ($(get(:prune_interval))) Iteration" : ""

            return pareto_str * cut_selection_str * pruning_str
        end

        return (; set!, get, get_setting_id)
    end
end
