using RCall, DataFrames 

include("RCall_Cut_conversion.jl")

function scalagon(cuts::Array{Cut})
    df = DataFrame(getfield.(cuts, :coefficients))
    df[!, :intercept] = getfield.(cuts, :intercept)
    pareto_indices = R"""
        library(rPref)
        # start_time = Sys.time()
        l = data.frame($df)
        # df_conversion_time = Sys.time() - start_time
        # print("Conversion time ist")
        # print(df_conversion_time)
        prefs = high("stock.3.")*  high("stock.2.")*  high("stock.1.")*  high("intercept")

        # start_time = Sys.time()
        # skyline = psel(l, prefs) 
        # psel_time = Sys.time() - start_time
        # print("PSEL time ist")
        # print(psel_time)

        # start_time = Sys.time()
        pareto_indices = psel.indices(l, prefs)
        # indices_time = Sys.time() - start_time
        # print("indices time ist")
        # print(indices_time)
        pareto_indices
    """
    return rcopy(Array{Number}, pareto_indices)
end