current_folder = pwd()

using Distributed

number_of_workers = nworkers() 
number_of_threads = Threads.nthreads()

# folder where data are stored
@eval @everywhere data_folder_path = "./" #$(ARGS[1])

# folder with fastwig tables to initialize the library
@eval @everywhere sl2cfoam_next_data_folder = "/home/frisus/Scrivania/sl2cfoam-next-dev/data_sl2cfoam" #$(ARGS[2])

printstyled("\nBlack-to-White hole amplitude parallelized on $(number_of_workers) worker(s)\n\n"; bold=true, color=:blue)

println("precompiling packages...")
@everywhere begin
    include("inc/pkgs.jl")
end
println("done\n")

println("precompiling source code...")
@everywhere begin
    include("parameters.jl")
    include("src/init.jl")
    include("src/utilities.jl")
    include("src/check.jl")
    include("src/vertex_functions.jl")
    include("src/generating_spins.jl")
end
println("done\n")

println("checking configurations to compute...")
check_on_preliminary_parameters(data_folder_path, sl2cfoam_next_data_folder, Dl_min, Dl_max, number_of_workers, number_of_threads)

@everywhere begin
    task_id = myid()
    number_of_tasks = nprocs()
    number_conf = size(angular_spins)[1]

    for user_conf in angular_spins
        check_configuration!(user_conf)
    end
end
println("done\n")

println("-------------------------------------------------------------------------\n")
printstyled("Starting computations\n\n"; bold=true, color=:blue)
println("-------------------------------------------------------------------------")

# Dataframes in which the computational times are stored 
computational_times = Array{Float64,2}(undef, Dl_max - Dl_min + 1, number_conf)
column_labels = String[]
counter_df = 0

for user_conf in angular_spins

    global counter_df += 1

    conf = init_config(user_conf, data_folder_path)
    @eval @everywhere conf = $conf

    printstyled("\n\nStarting with configuration:\nj0=$(conf.j0), jpm=$(conf.jpm) ...\n\n"; bold=true, color=:bold)
    sleep(1)


    #####################################################################################################################################
    # COMPUTING SPINS COMBINATIONS
    #####################################################################################################################################

    printstyled("\nComputing all possible spins combinations\n",
        "for angular spins in [$(conf.j0_float - conf.K0), $(conf.j0_float + conf.K0)] ",
        "and radial spins in [$(conf.jpm_float - conf.Kpm), $(conf.jpm_float + conf.Kpm)] ...\n\n"; bold=true, color=:cyan)
    generating_spins(conf.j0, conf.K0, conf.jpm, conf.Kpm, conf.base_folder, conf.spinfoam_folder, immirzi, verbosity_flux)

    #####################################################################################################################################
    # COMPUTING VERTICES
    #####################################################################################################################################


end


current_date = now()
pre_final_data_path = "$(current_folder)/data"
final_data_path = "$(pre_final_data_path)/immirzi_$(immirzi)/alpha_$(alpha)"
mkpath(final_data_path)