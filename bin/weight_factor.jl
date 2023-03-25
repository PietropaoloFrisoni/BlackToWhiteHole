current_folder = pwd()

using Distributed

@eval @everywhere number_of_workers = nworkers()
@eval @everywhere number_of_threads = Threads.nthreads()

# folder where data are stored
@eval @everywhere data_folder_path = $(ARGS[1])

# folder with fastwigxj tables to initialize the library
@eval @everywhere sl2cfoam_next_data_folder = $(ARGS[2])

printstyled("\nBlack-to-White hole amplitude parallelized on $(number_of_workers) worker(s)\n\n"; bold=true, color=:blue)

println("precompiling packages...")
@everywhere begin
    include("../inc/pkgs.jl")
end
println("done\n")

println("precompiling source code...")
@everywhere begin
    include("../parameters.jl")
    include("../src/init.jl")
    include("../src/utilities.jl")
    include("../src/check.jl")
    include("../src/amplitude.jl")
end
println("done\n")

println("checking configurations to compute...")
CheckPreliminaryParameters(data_folder_path, sl2cfoam_next_data_folder, Dl_min, Dl_max, number_of_workers, number_of_threads)

@everywhere begin
    task_id = myid()
    number_of_tasks = nprocs()
    number_conf = size(angular_spins)[1]

    for user_conf in angular_spins
        CheckConfiguration!(user_conf)
    end
end
println("done\n")

println("-------------------------------------------------------------------------\n")
printstyled("Starting computations\n\n"; bold=true, color=:blue)
println("-------------------------------------------------------------------------")

for user_conf in angular_spins

    local conf = InitConfig(user_conf, data_folder_path)
    @eval @everywhere conf = $conf

    printstyled("\n\nStarting with configuration:\nj0=$(conf.j0), jpm=$(conf.jpm) ...\n\n"; bold=true, color=:bold)
    sleep(1)

    #####################################################################################################################################
    ### COMPUTING THE WEIGHT FACTOR
    #####################################################################################################################################

    printstyled("\nComputing the weight factor...\n"; bold=true, color=:blue)

    @everywhere begin
        @load "$(conf.base_folder)/spins_configurations.jld2" spins_configurations
        @load "$(conf.base_folder)/spins_map.jld2" spins_map
        @load "$(conf.base_folder)/intertwiners_range.jld2" intertwiners_range
    end

    m = sqrt(conf.j0_float * immirzi)
    T_range = LinRange(0, 4 * pi * m / immirzi, T_sampling_parameter)

    @time @sync @distributed for current_angular_spins_comb in eachindex(spins_map)

        total_radial_spins_combinations = spins_map[current_angular_spins_comb]
        upper_bound = sum(spins_map[1:current_angular_spins_comb])
        lower_bound = upper_bound - total_radial_spins_combinations + 1

        j1 = twice(spins_configurations[lower_bound][1]) / 2
        j2 = twice(spins_configurations[lower_bound][2]) / 2
        j3 = twice(spins_configurations[lower_bound][3]) / 2
        j4 = twice(spins_configurations[lower_bound][4]) / 2

        path_weight_factor = "$(conf.spinfoam_folder)/j1_$(j1)_j2_$(j2)_j3_$(j3)_j4_$(j4)/immirzi_$(immirzi)/alpha_$(alpha)"
        mkpath("$(path_weight_factor)")

        total_elements = total_radial_spins_combinations^2

        weight_factor = Array{ComplexF64,2}(undef, total_elements, T_sampling_parameter)

        WeightFactor!(weight_factor, alpha, conf.j0_float, conf.jpm_float, m, T_range, immirzi, spins_configurations,
            lower_bound, upper_bound, j1, j2, j3, j4)

        @save "$(path_weight_factor)/weight_factor_T_$(T_sampling_parameter).jld2" weight_factor

    end

end

# release workers
if (number_of_workers > 1)
    for i in workers()
        rmprocs(i)
    end
end

printstyled("\nCompleted\n\n"; bold=true, color=:blue)