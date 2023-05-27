current_folder = pwd()

@eval @everywhere number_of_workers = nworkers()
@eval @everywhere number_of_threads = Threads.nthreads()

# folder where data are stored
@eval @everywhere data_folder_path = $(ARGS[1])

# folder with fastwigxj tables to initialize the library
@eval @everywhere sl2cfoam_next_data_folder = $(ARGS[2])

printstyled("\nBlack-to-White hole lifetime computation parallelized on $(number_of_workers) worker(s)\n\n"; bold=true, color=:blue)

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
    number_conf = size(angular_spins,1)

    for user_conf in angular_spins
        CheckConfiguration!(user_conf)
    end
end
verbosity_flux && println("done\n")

println("-------------------------------------------------------------------------\n")
printstyled("Starting computations\n\n"; bold=true, color=:blue)
println("-------------------------------------------------------------------------")

for user_conf in angular_spins

    local conf = InitConfig(user_conf, data_folder_path)
    @eval @everywhere conf = $conf

    printstyled("\n\nStarting with configuration:\nj0=$(conf.j0_float), jpm=$(conf.jpm_float) ...\n\n"; bold=true, color=:bold)
    sleep(1)

    

end

ReleaseWorkers(number_of_workers)
printstyled("\nCompleted\n\n"; bold=true, color=:blue)