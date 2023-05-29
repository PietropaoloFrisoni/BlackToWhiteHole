current_folder = pwd()

@eval @everywhere number_of_workers = nworkers()
@eval @everywhere number_of_threads = Threads.nthreads()

# folder where data are stored
@eval @everywhere data_folder_path = $(ARGS[1])

# folder with fastwigxj tables to initialize the library
@eval @everywhere sl2cfoam_next_data_folder = $(ARGS[2])

printstyled("\nBlack-to-White hole amplitude assembling parallelized on $(number_of_workers) worker(s)\n\n"; bold=true, color=:blue)

@everywhere include("../parameters.jl")

verbosity_flux && println("precompiling packages...")
@everywhere begin
    include("../inc/pkgs.jl")
end
verbosity_flux && println("done\n")

verbosity_flux && println("precompiling source code...")
@everywhere begin
    include("../parameters.jl")
    include("../src/init.jl")
    include("../src/utilities.jl")
    include("../src/check.jl")
    include("../src/amplitude.jl")
end
verbosity_flux && println("done\n")

verbosity_flux && println("checking configurations to compute...")
CheckPreliminaryParameters(data_folder_path, sl2cfoam_next_data_folder, Dl_min, Dl_max, 1, 0)

@everywhere begin
    task_id = myid()
    number_of_tasks = nprocs()
    number_conf = size(angular_spins, 1)

    for user_conf in angular_spins
        CheckConfiguration!(user_conf)
    end
end
verbosity_flux && println("done\n")

number_of_workers > T_sampling_parameter && warn("Parallelization not fully exploited: T_sampling_parameter is too low")

println("-------------------------------------------------------------------------\n")
printstyled("Starting computations\n\n"; bold=true, color=:blue)
println("-------------------------------------------------------------------------")

for user_conf in angular_spins

    local conf = InitConfig(user_conf, data_folder_path)
    @eval @everywhere conf = $conf

    printstyled("\n\nStarting with configuration:\nj0=$(conf.j0_float), jpm=$(conf.jpm_float) ...\n\n"; bold=true, color=:bold)
    sleep(1)

    #####################################################################################################################################
    ### ASSEMBLING THE AMPLITUDE
    #####################################################################################################################################

    printstyled("\nAssembling the amplitude...\n"; bold=true, color=:blue)

    @everywhere begin
        @load "$(conf.base_folder)/spins_configurations.jld2" spins_configurations
        @load "$(conf.base_folder)/spins_map.jld2" spins_map
        @load "$(conf.base_folder)/intertwiners_range.jld2" intertwiners_range
    end

    T_range = LinRange(0, 4 * pi * conf.m / immirzi, T_sampling_parameter)
    Dl_range = Dl_max - Dl_min + 1

    amplitude = SharedArray{ComplexF64}(T_sampling_parameter, Dl_range)
    amplitude[:] .= 0.0 + 0.0im

    local column_labels = String[]

    for Dl = Dl_min:Dl_max

        printstyled("\nCurrent Dl=$(Dl)...\n"; bold=true, color=:magenta)
        push!(column_labels, "Dl_$(Dl)")

        @time @sync @distributed for T = 1:T_sampling_parameter

            for current_angular_spins_comb in eachindex(spins_map)

                total_radial_spins_combinations = spins_map[current_angular_spins_comb]
                upper_bound = sum(spins_map[1:current_angular_spins_comb])
                lower_bound = upper_bound - total_radial_spins_combinations + 1

                j1 = twice(spins_configurations[lower_bound][1]) / 2
                j2 = twice(spins_configurations[lower_bound][2]) / 2
                j3 = twice(spins_configurations[lower_bound][3]) / 2
                j4 = twice(spins_configurations[lower_bound][4]) / 2

                path_contracted_spinfoam = "$(conf.spinfoam_folder)/j1_$(j1)_j2_$(j2)_j3_$(j3)_j4_$(j4)/immirzi_$(immirzi)/Dl_$(Dl)"
                path_weight_factor = "$(conf.spinfoam_folder)/j1_$(j1)_j2_$(j2)_j3_$(j3)_j4_$(j4)/immirzi_$(immirzi)/alpha_$(alpha)"

                @load "$(path_contracted_spinfoam)/contracted_spinfoam.jld2" contracted_spinfoam
                @load "$(path_weight_factor)/weight_factor_T_$(T_sampling_parameter).jld2" weight_factor

                amplitude[T, Dl+1] += dot(weight_factor[:, T], contracted_spinfoam[:])

            end

        end

    end

    amplitude_abs_sq = zeros(T_sampling_parameter, Dl_range)
    amplitude_abs_sq_integrated = zeros(Dl_range)
    amplitude_abs_sq_T_integrated = zeros(Dl_range)

    [[amplitude_abs_sq[T, Dl+1] = abs(amplitude[T, Dl+1])^2 * GlobalFactor(conf.j0_float, conf.jpm_float, alpha) for T = 1:T_sampling_parameter] for Dl = Dl_min:Dl_max]
    sum_check(amplitude_abs_sq) && error("NaN or Inf in amplitude absolute squared")

    AmplitudeIntegration!(amplitude_abs_sq_integrated, amplitude_abs_sq_T_integrated, amplitude_abs_sq, T_range, T_sampling_parameter)

    crossing_times = [amplitude_abs_sq_T_integrated[Dl_index] / amplitude_abs_sq_integrated[Dl_index] for Dl_index in 1:Dl_range]
    sum_check(crossing_times) && error("NaN or Inf in crossing times")

    amplitude_abs_sq_df = DataFrame(amplitude_abs_sq, column_labels)
    crossing_times_df = DataFrame(transpose(crossing_times), column_labels)

    base_folder_alpha = "$(conf.base_folder)/immirzi_$(immirzi)/alpha_$(alpha)/Dl_range_$(Dl_range)"
    mkpath(base_folder_alpha)
    CSV.write("$(base_folder_alpha)/amplitude_abs_sq_T_$(T_sampling_parameter).csv", amplitude_abs_sq_df)
    CSV.write("$(base_folder_alpha)/crossing_time_$(T_sampling_parameter).csv", crossing_times_df)

end

ReleaseWorkers(number_of_workers)
printstyled("\nCompleted\n\n"; bold=true, color=:blue)