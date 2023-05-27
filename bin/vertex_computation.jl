current_folder = pwd()

using Distributed

@eval @everywhere number_of_workers = nworkers()
@eval @everywhere number_of_threads = Threads.nthreads()

# folder where data are stored
@eval @everywhere data_folder_path = $(ARGS[1])

# folder with fastwigxj tables to initialize the library
@eval @everywhere sl2cfoam_next_data_folder = $(ARGS[2])

printstyled("\nBlack-to-White hole vertex computations parallelized on $(number_of_workers) worker(s)\n\n"; bold=true, color=:blue)

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
    include("../src/vertex_functions.jl")
    include("../src/generating_spins.jl")
end
verbosity_flux && println("done\n")

verbosity_flux && println("checking configurations to compute...")
CheckPreliminaryParameters(data_folder_path, sl2cfoam_next_data_folder, Dl_min, Dl_max, number_of_workers, number_of_threads)

@everywhere begin
    task_id = myid()
    number_of_tasks = nprocs()
    number_conf = size(angular_spins,1)

    for user_conf in angular_spins
        CheckConfiguration!(user_conf)
    end
end
verbosity_flux && println("done\n")

verbosity_flux && log("Initializing sl2cfoam-next on each worker...")
@everywhere InitSL2Cfoam(immirzi, sl2cfoam_next_data_folder, number_of_threads, verbosity_flux)
verbosity_flux && println("done\n")

current_date = now()
data_folder_path = "$(data_folder_path)"
comp_times_data_path = "$(data_folder_path)/data/computational_times/run_started_on:$(current_date)"
mkpath(comp_times_data_path)

println("-------------------------------------------------------------------------\n")
printstyled("Starting computations\n\n"; bold=true, color=:blue)
println("-------------------------------------------------------------------------")

# Dataframe in which the computational times are stored 
computational_times = Array{Float64,2}(undef, Dl_max - Dl_min + 1, number_conf)
column_labels = String[]
counter_df = 0

for user_conf in angular_spins

    global counter_df += 1

    conf = InitConfig(user_conf, data_folder_path)
    @eval @everywhere conf = $conf

    printstyled("\n\nStarting with configuration:\nj0=$(conf.j0), jpm=$(conf.jpm) ...\n\n"; bold=true, color=:bold)
    sleep(1)

    #####################################################################################################################################
    ### COMPUTING SPINS COMBINATIONS
    #####################################################################################################################################

    printstyled("\nComputing spins combinations\n",
        "for angular spins in [$(conf.j0_float - conf.K0), $(conf.j0_float + conf.K0)] ",
        "and radial spins in [$(conf.jpm_float - conf.Kpm), $(conf.jpm_float + conf.Kpm)] ...\n\n"; bold=true, color=:cyan)
    GeneratingSpins(conf.j0, conf.K0, conf.jpm, conf.Kpm, conf.base_folder, conf.spinfoam_folder, immirzi, verbosity_flux)

    #####################################################################################################################################
    ### COMPUTING VERTICES
    #####################################################################################################################################

    @load "$(conf.base_folder)/spins_configurations.jld2" spins_configurations

    printstyled("\nComputing all vertices...\n"; bold=true, color=:blue)

    for Dl = Dl_min:Dl_max

        printstyled("\nCurrent Dl=$(Dl)...\n"; bold=true, color=:magenta)

        seconds_required = @elapsed VertexDistributeSingleMachine(spins_configurations, Dl, (true, true))
        println("done. Time required: $(seconds_required) seconds\n")
        computational_times[Dl+1, counter_df] = seconds_required

    end

    push!(column_labels, "j0 = $(conf.j0_float), jpm = $(conf.jpm_float)")

end

# store computational times
df = DataFrame(computational_times, column_labels)
CSV.write("$(comp_times_data_path)/immirzi_$(immirzi)_workers_$(number_of_workers)_threads_$(number_of_threads).csv", df)

printstyled("\nCompleted\n\n"; bold=true, color=:blue)