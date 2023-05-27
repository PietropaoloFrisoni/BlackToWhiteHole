current_folder = pwd()

using Distributed

@eval @everywhere number_of_workers = nworkers()
@eval @everywhere number_of_threads = Threads.nthreads()

# folder where data are stored
@eval @everywhere data_folder_path = $(ARGS[1])

# folder with fastwigxj tables to initialize the library
@eval @everywhere sl2cfoam_next_data_folder = $(ARGS[2])

printstyled("\nBlack-to-White hole vertex contraction parallelized on $(number_of_workers) worker(s)\n\n"; bold=true, color=:blue)

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
    include("../src/vertex_functions.jl")
    include("../src/generating_spins.jl")
    include("../src/wigner_matrix.jl")
end
println("done\n")

println("checking configurations to compute...")
CheckPreliminaryParameters(data_folder_path, sl2cfoam_next_data_folder, Dl_min, Dl_max, number_of_workers, number_of_threads)

@everywhere begin
    task_id = myid()
    number_of_tasks = nprocs()
    number_conf = size(angular_spins,1)

    for user_conf in angular_spins
        CheckConfiguration!(user_conf)
    end
end
println("done\n")

println("Initializing sl2cfoam-next on each worker...")
@everywhere InitSL2Cfoam(immirzi, sl2cfoam_next_data_folder, number_of_threads, verbosity_flux)
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
    ### CONTRACTING VERTICES
    #####################################################################################################################################

    printstyled("\nContracting vertices...\n"; bold=true, color=:blue)

    @everywhere begin
        @load "$(conf.base_folder)/spins_configurations.jld2" spins_configurations
        @load "$(conf.base_folder)/spins_map.jld2" spins_map
        @load "$(conf.base_folder)/intertwiners_range.jld2" intertwiners_range
    end

    for Dl = Dl_min:Dl_max

        printstyled("\nCurrent Dl=$(Dl)...\n"; bold=true, color=:magenta)

        @time @sync @distributed for current_angular_spins_comb in eachindex(spins_map)

            total_radial_spins_combinations = spins_map[current_angular_spins_comb]
            upper_bound = sum(spins_map[1:current_angular_spins_comb])
            lower_bound = upper_bound - total_radial_spins_combinations + 1

            j1 = twice(spins_configurations[lower_bound][1]) / 2
            j2 = twice(spins_configurations[lower_bound][2]) / 2
            j3 = twice(spins_configurations[lower_bound][3]) / 2
            j4 = twice(spins_configurations[lower_bound][4]) / 2

            path_contracted_spinfoam = "$(conf.spinfoam_folder)/j1_$(j1)_j2_$(j2)_j3_$(j3)_j4_$(j4)/immirzi_$(immirzi)/Dl_$(Dl)"
            mkpath(path_contracted_spinfoam)

            i1_range = intertwiners_range[lower_bound][1]
            total_elements = total_radial_spins_combinations^2

            coherent_matrix_up = Array{ComplexF64,2}(undef, i1_range, total_radial_spins_combinations)
            coherent_matrix_down = Array{ComplexF64,2}(undef, i1_range, total_radial_spins_combinations)
            contracted_spinfoam = Vector{ComplexF64}(undef, total_elements)

            SpinfoamContractUp!(coherent_matrix_up, lower_bound, upper_bound, i1_range, spins_configurations, Dl)

            # Since there is complete symmetry between top and bottom, we can perform the contraction just once
            # SpinfoamContractDown!(coherent_matrix_down, lower_bound, upper_bound, i1_range, spins_configurations, Dl)
            coherent_matrix_down[:] .= coherent_matrix_up[:]

            SpinfoamFinalContraction!(contracted_spinfoam, coherent_matrix_up, coherent_matrix_down, i1_range, total_radial_spins_combinations)
            sum_check(contracted_spinfoam) && error("NaN or Inf in contracted spinfoam")

            @save "$(path_contracted_spinfoam)/contracted_spinfoam.jld2" contracted_spinfoam

        end

    end

end

printstyled("\nCompleted\n\n"; bold=true, color=:blue)