# Execute the upstream pixel fitter and decoder without modifying either package.
using SeqFISH_ADCG, SeqFISHSyndromeDecoding, CSV, DataFrames, GLPK, JSON3
using DelimitedFiles, LinearAlgebra, Random, Pkg, SHA, TOML

function main(directory, destination)
    Random.seed!(2026110106)
    BLAS.set_num_threads(1)
    input_bytes = read(joinpath(directory, "inputs.json"))
    input = JSON3.read(String(copy(input_bytes)))
    cb = CSV.read(joinpath(directory, "codebook.csv"), DataFrame)
    H = Int.(readdlm(joinpath(directory, "parity_check.csv"), ',', Int))
    @assert all(H * Matrix(cb[:, 2:end])' .% 2 .== 0)
    params = DecodeParams()
    set_zeros_probed(params, false)
    set_n_allowed_drops(params, 0)
    set_xy_search_radius(params, 4)
    set_z_search_radius(params, 0)
    set_lat_var_cost_coeff(params, 7)
    set_lw_var_cost_coeff(params, 3.2)
    set_s_var_cost_coeff(params, 0)
    set_z_var_cost_coeff(params, 0)
    set_free_dot_cost(params, 1)
    outputs = []
    for stack in input.stacks
        points = DataFrame(x=Float64[], y=Float64[], s=Float64[], w=Float64[], z=Float64[], block=UInt8[], pseudocolor=UInt8[])
        for (bit, image) in enumerate(stack.images)
            tile = Float64.(permutedims(reduce(hcat, image)))
            records = fit_tile((tile, 1., 1.4, Float64(stack.background), 1.e8, 1., Float64(stack.min_weight), Int(stack.max_iters), 20, "ADCG"))
            fitted = records.last_iteration
            for row in eachrow(fitted)
                push!(points, (row.x-1, row.y-1, row.s, row.w, 0., UInt8(bit), 0x01))
            end
            checkpoint = Dict("stack"=>String(stack.key), "completed_bits"=>bit,
                              "fitted_points"=>[Dict(string(k)=>row[k] for k in propertynames(points)) for row in eachrow(points)],
                              "completed_stack_outputs"=>outputs)
            write(destination*".checkpoint.json", JSON3.write(checkpoint))
            println(stack.key, " bit ", bit, ": ", nrow(fitted), " fitted points")
            flush(stdout)
        end
        SeqFISHSyndromeDecoding.sort_readouts!(points)
        saved_points = [Dict(string(k)=>row[k] for k in propertynames(points)) for row in eachrow(points)]
        decoded = nrow(points) == 0 ? nothing : decode_syndromes!(points, cb, H, params; optimizer=GLPK.Optimizer)
        genes = isnothing(decoded) ? String[] : String.(decoded.gene)
        result = Dict("key"=>String(stack.key), "mode"=>String(stack.mode), "case"=>stack.case,
                      "fitted_points"=>saved_points, "decoded_genes"=>genes,
                      "decoder_empty"=>isnothing(decoded))
        result["saturated_bits"] = [bit for bit in 1:16 if sum(points.block .== bit) >= stack.max_iters]
        if stack.mode != "measured"
            first, second = String.(stack.first), String.(stack.alternative)
            score = sum(g in first for g in genes)-sum(g in second for g in genes)
            result["claim_score"] = score
            result["forced_decision"] = score < 0 ? 1 : 0
            if stack.mode == "alias"
                result["paired_errors"] = [result["forced_decision"] != 0, result["forced_decision"] != 1]
            end
            if stack.mode != "alias"
                truth = stack.mode == "registered_first" ? first : second
                result["exact_pair_recovered"] = sort(genes) == sort(truth)
            end
        end
        push!(outputs, result)
        println(stack.key, ": ", nrow(points), " dots; ", join(genes, ", "))
        flush(stdout)
    end
    deps = [Dict("name"=>info.name, "version"=>string(info.version), "tree_hash"=>string(info.tree_hash))
            for info in values(Pkg.dependencies()) if info.name in ["SeqFISH_ADCG", "SeqFISHSyndromeDecoding", "GLPK"]]
    manifest = TOML.parsefile(joinpath(dirname(Base.active_project()), "Manifest.toml"))["deps"]
    commits = Dict(name=>manifest[name][1]["repo-rev"] for name in ["SeqFISH_ADCG", "SeqFISHSyndromeDecoding"])
    controls = [r["exact_pair_recovered"] for r in outputs if startswith(r["mode"], "registered")]
    result = Dict("protocol_commit"=>String(input.protocol_commit), "julia_version"=>string(VERSION),
                  "packages"=>deps, "outputs"=>outputs, "registered_control_accuracy"=>sum(controls)/length(controls),
                  "inputs_sha256"=>bytes2hex(sha256(input_bytes)), "package_commits"=>commits,
                  "gates"=>Dict("registered_control"=>sum(controls)/length(controls)>=.8,
                                "paired_alias_error"=>all(sum(r["paired_errors"]) == 1 for r in outputs if r["mode"] == "alias")),
                  "independent_biological_validation"=>false, "poisson_noise_evaluated"=>false)
    write(destination, JSON3.write(result))
end

main(ARGS[1], ARGS[2])
