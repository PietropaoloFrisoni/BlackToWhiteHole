current_folder = pwd()

@eval number_of_workers = 1
@eval number_of_threads = 0

# folder where data are stored
@eval data_folder_path = $(ARGS[1])

# folder with fastwigxj tables to initialize the library
@eval sl2cfoam_next_data_folder = $(ARGS[2])

# number of points to plot amplitude as function of T
@eval T_sampling_parameter = parse(Int, ARGS[3])

printstyled("\nBlack-to-White hole amplitude\n\n"; bold=true, color=:blue)

println("precompiling packages...")
include("../inc/pkgs.jl")
println("done\n")

println("precompiling source code...")
include("../parameters.jl")
include("../src/init.jl")
include("../src/utilities.jl")
include("../src/check.jl")
include("../src/amplitude.jl")
println("done\n")

println("checking configurations to compute...")
CheckPreliminaryParameters(data_folder_path, sl2cfoam_next_data_folder, Dl_min, Dl_max, 1, 0)

number_conf = size(angular_spins)[1]

for user_conf in angular_spins
    CheckConfiguration!(user_conf)
end
println("done\n")

println("-------------------------------------------------------------------------\n")
printstyled("Starting computations\n\n"; bold=true, color=:blue)
println("-------------------------------------------------------------------------")

for user_conf in angular_spins

    conf = InitConfig(user_conf, data_folder_path)
    @eval conf = $conf

    printstyled("\n\nStarting with configuration:\nj0=$(conf.j0), jpm=$(conf.jpm) ...\n\n"; bold=true, color=:bold)
    sleep(1)

    number_of_T_points = T_sampling_parameter

    #####################################################################################################################################
    ### ASSEMBLING THE AMPLITUDE
    #####################################################################################################################################

    printstyled("\nAssembling the amplitude...\n"; bold=true, color=:blue)

    @load "$(conf.base_folder)/spins_configurations.jld2" spins_configurations
    @load "$(conf.base_folder)/spins_map.jld2" spins_map
    @load "$(conf.base_folder)/intertwiners_range.jld2" intertwiners_range

    m = sqrt(conf.j0_float * immirzi)
    T_range = LinRange(0, 4 * pi * m / immirzi, T_sampling_parameter)

    amplitude = Array{ComplexF64,2}(undef, number_of_T_points, Dl_max - Dl_min + 1)
    amplitude_abs_sq = zeros(number_of_T_points, Dl_max - Dl_min + 1) #Array{Float64,2}(undef, number_of_T_points, Dl_max - Dl_min + 1)
    amplitude_abs_sq_integrated = zeros(Dl_max - Dl_min + 1) #Array{Float64,1}(undef, Dl_max - Dl_min + 1)
    #amplitude_abs_sq_integrated[:] .= 0.0
    amplitude[:] .= 0.0 + 0.0 * im
    #amplitude_abs_sq[:] .= 0.0

    column_labels = String[]

    for Dl = Dl_min:Dl_max

        printstyled("\nCurrent Dl=$(Dl)...\n"; bold=true, color=:magenta)

        push!(column_labels, "Dl_$(Dl)")

        @time for T = 1:number_of_T_points

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
                @load "$(path_weight_factor)/weight_factor_T_$(number_of_T_points).jld2" weight_factor

                amplitude[T, Dl+1] += dot(weight_factor[:, T], contracted_spinfoam[:])

            end

        end

    end

    for Dl = Dl_min:Dl_max
        for T = 1:number_of_T_points
            amplitude_abs_sq[T, Dl+1] = abs(amplitude[T, Dl+1])^2
        end
    end

    amplitude_abs_sq_df = DataFrame(amplitude_abs_sq, column_labels)

    


end