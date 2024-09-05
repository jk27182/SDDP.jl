mutable struct Settings
    use_pareto_cut_logic::Bool
    log_level::Int
    use_cut_selection::Bool
    use_pruning::Bool
    prune_interval::Int
    debug_mode::Bool
    problem_name::String
    cut_type::Union{Enum, Nothing} 

    # use a closure to store the settings and prevent public acces
    function Settings(;
        use_pareto_cut_logic::Bool,
        log_level::Int,
        use_cut_selection::Bool,
        use_pruning::Bool,
        prune_interval::Int,
        debug_mode::Bool,
        problem_name::String="",
        cut_type::Union{Enum, Nothing}=nothing,
    )
        settings = new(
            use_pareto_cut_logic,
            log_level,
            use_cut_selection,
            use_pruning,
            prune_interval,
            debug_mode,
            problem_name,
            cut_type,
        )
        function set!(;
            use_pareto_cut_logic::Bool,
            log_level::Int,
            use_cut_selection::Bool,
            use_pruning::Bool,
            prune_interval::Int,
            debug_mode::Bool,
            problem_name::String,
            cut_type::Enum,   
        )
            settings.use_pareto_cut_logic = use_pareto_cut_logic
            settings.log_level = log_level
            settings.use_cut_selection = use_cut_selection
            settings.use_pruning = use_pruning
            settings.prune_interval = prune_interval
            settings.debug_mode = debug_mode
            settings.problem_name = problem_name
            settings.cut_type = cut_type    
            return
        end
        set!(field::Symbol, value) = setfield!(settings, field, value)
        set!(field::String, value) = setfield!(settings, Symbol(field), value)

        get() = settings
        get(x::Symbol) = getfield(settings, x)
        get(x::String) = getfield(settings, Symbol(x))

        function get_setting_id()
            cut_str = get(:cut_type) == SINGLE_CUT ? "SingleCut_" : "MultiCut_"
            pareto_str = get(:use_pareto_cut_logic) ? "ParetoCuts_" : ""
            cut_selection_str = get(:use_cut_selection) ? "StandardCutSelection_" : ""
            pruning_str = get(:use_pruning) ?  "Pruning$(get(:prune_interval))" : ""

            id_str = pareto_str * cut_selection_str * pruning_str
            id_str = isempty(id_str) ? "DefaultSDDP" : id_str
            return id_str * cut_str
        end

        return (; set!, get, get_setting_id)
    end
end
