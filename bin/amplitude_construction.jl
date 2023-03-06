current_folder = pwd()

using Distributed

@eval @everywhere number_of_workers = nworkers()
@eval @everywhere number_of_threads = Threads.nthreads()

# folder where data are stored
@eval @everywhere data_folder_path = $(ARGS[1])

# folder with fastwigxj tables to initialize the library
@eval @everywhere sl2cfoam_next_data_folder = $(ARGS[2])

# number of points to plot amplitude as function of T
@eval @everywhere T_sampling_parameter = parse(Int, ARGS[3])

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
    include("../src/vertex_functions.jl")
    include("../src/amplitude.jl")
    include("../src/generating_spins.jl")
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

println("Initializing sl2cfoam-next on each worker...")
@everywhere InitSL2Cfoam(immirzi, sl2cfoam_next_data_folder, number_of_threads, verbosity_flux)
println("done\n")

current_date = now()
comp_times_data_path = "$(data_folder_path)/data/computational_times/vertex_computations/run_started_on:$(current_date)"
mkpath(comp_times_data_path)

println("-------------------------------------------------------------------------\n")
printstyled("Starting computations\n\n"; bold=true, color=:blue)
println("-------------------------------------------------------------------------")

for user_conf in angular_spins

    conf = InitConfig(user_conf, data_folder_path)
    @eval @everywhere conf = $conf

    printstyled("\n\nStarting with configuration:\nj0=$(conf.j0), jpm=$(conf.jpm) ...\n\n"; bold=true, color=:bold)
    sleep(1)

    #####################################################################################################################################
    ### ASSEMBLING BLACK-TO-WHITE HOLE AMPLITUDE
    #####################################################################################################################################

    printstyled("\nAssembling B-W amplitude...\n"; bold=true, color=:blue)

    @everywhere begin
        @load "$(conf.base_folder)/spins_configurations.jld2" spins_configurations
        @load "$(conf.base_folder)/spins_map.jld2" spins_map
        @load "$(conf.base_folder)/intertwiners_range.jld2" intertwiners_range
    end

    m = sqrt(conf.j0_float * immirzi)
    T_range = LinRange(0, 4 * pi * m / immirzi, T_sampling_parameter)
    number_of_T_points = size(T_range)[1]

    for Dl = Dl_min:Dl_max

        printstyled("\nCurrent Dl=$(Dl)...\n"; bold=true, color=:magenta)

        # Dataframes in which the amplitudes will be stored at the end       
        BW_amplitudes = Array{ComplexF64,2}(undef, number_of_T_points + 1, number_conf)
        BW_amplitudes_abs_sq = Array{Float64,2}(undef, number_of_T_points + 1, number_conf)
        BW_amplitudes_abs_sq_integrated = Array{Float64,2}(undef, 1, number_conf)

        @time @sync @distributed for current_angular_spins_comb in eachindex(spins_map)

            total_radial_spins_combinations = spins_map[current_angular_spins_comb]
            upper_bound = sum(spins_map[1:current_angular_spins_comb])
            lower_bound = upper_bound - total_radial_spins_combinations + 1

            j1 = twice(spins_configurations[lower_bound][1]) / 2
            j2 = twice(spins_configurations[lower_bound][2]) / 2
            j3 = twice(spins_configurations[lower_bound][3]) / 2
            j4 = twice(spins_configurations[lower_bound][4]) / 2

            path_contracted_spinfoam = "$(conf.spinfoam_folder)/j1_$(j1)_j2_$(j2)_j3_$(j3)_j4_$(j4)/immirzi_$(immirzi)/Dl_$(Dl)"
            path_amplitude = "$(path_contracted_spinfoam)/alpha_$(alpha)"
            mkpath("$(path_amplitude)")
            @load "$(path_contracted_spinfoam)/contracted_spinfoam.jld2" contracted_spinfoam

            i1_range = intertwiners_range[lower_bound][1]
            total_elements = total_radial_spins_combinations^2

            amplitude = Array{ComplexF64,1}(undef, number_of_T_points)

            get_T_points!(amplitude, contracted_spinfoam, alpha, conf.j0_float, conf.jpm_float, m, T_range, immirzi, spins_configurations,
                lower_bound, upper_bound, path_contracted_spinfoam, j1, j2, j3, j4)

        end

        printstyled("\nSumming over angular spins...\n"; bold=true, color=:cyan)

        final_amplitude = Array{ComplexF64,1}(undef, number_of_T_points)
        final_amplitude_abs_sq = Array{Float64,1}(undef, number_of_T_points)

        for T_index in eachindex(T_range)
            final_amplitude[T_index] = 0.0 + 0.0 * im
            final_amplitude_abs_sq[T_index] = 0.0
        end

        @time @inbounds for current_angular_spins_comb in eachindex(spins_map)

            total_radial_spins_combinations = spins_map[current_angular_spins_comb]
            upper_bound = sum(spins_map[1:current_angular_spins_comb])
            lower_bound = upper_bound - total_radial_spins_combinations + 1

            j1 = twice(spins_configurations[lower_bound][1]) / 2
            j2 = twice(spins_configurations[lower_bound][2]) / 2
            j3 = twice(spins_configurations[lower_bound][3]) / 2
            j4 = twice(spins_configurations[lower_bound][4]) / 2

            path_contracted_spinfoam = "$(conf.spinfoam_folder)/j1_$(j1)_j2_$(j2)_j3_$(j3)_j4_$(j4)/immirzi_$(immirzi)/Dl_$(Dl)"
            path_amplitude = "$(path_contracted_spinfoam)/alpha_$(alpha)"
            @load "$(path_amplitude)/amplitude_T_$(number_of_T_points).jld2" amplitude

            for T_index in eachindex(T_range)

                if (isnan(amplitude[T_index]))
                    println("NaN for con with j1_$(j1)_j2_$(j2)_j3_$(j3)_j4_$(j4) at index $(T_index)")
                end

                final_amplitude[T_index] += amplitude[T_index]

            end

        end

    end

end

# release workers
if (number_of_workers > 1)
    for i in workers()
        rmprocs(i)
    end
end

printstyled("\nCompleted\n\n"; bold=true, color=:blue)