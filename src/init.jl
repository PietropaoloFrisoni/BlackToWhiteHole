mutable struct Configuration

    j0::Half{Int8}
    j0_float::Float64
    K0::Half{Int8}
    jpm::Half{Int8}
    jpm_float::Float64
    Kpm::Half{Int8}
    m::Float64

    contracted_spinfoam_found::Bool

    base_folder::String
    spinfoam_folder::String

end

function InitConfig(user_conf, data_folder_path::String, contracted_spinfoam_found::Bool=false)

    j0_float = user_conf[1]
    jpm_float = round(Int64, 2 * (j0_float / (sqrt(6)))) / 2
    j0 = half(2 * j0_float)
    jpm = half(2 * jpm_float)

    K0 = half(1)
    Kpm = half(1)

    # where spins_conf, spins_map and intertwiner_range are stored
    base_folder = "$(data_folder_path)/data/amplitude_data/j0=$(j0_float)_jpm=$(jpm_float)/K0_$(twice(K0)/2)_Kpm_$(twice(Kpm)/2)"

    # subfolder with all angular spins combinations
    spinfoam_folder = "$(base_folder)/angular_spins"

    # This estimate disregards the term exp(-T/2m)
    m = sqrt(2 * immirzi * j0_float)

    Configuration(j0, j0_float, K0, jpm, jpm_float, Kpm, m, contracted_spinfoam_found, base_folder, spinfoam_folder)

end

function InitSL2Cfoam(immirzi, sl2cfoam_next_data_folder::String, number_of_threads::Int64, verbosity_flux::Bool; shell_parallelization::Bool=false)

    isMPI = @ccall SL2Cfoam.clib.sl2cfoam_is_MPI()::Bool
    isMPI && error("MPI version not allowed")

    conf_sl2cfoam_next = SL2Cfoam.Config(VerbosityOff, HighAccuracy, 100, 0)
    SL2Cfoam.cinit(sl2cfoam_next_data_folder, immirzi, conf_sl2cfoam_next)

    # enable C library automatic parallelization
    number_of_threads > 1 && (shell_parallelization = true)
    SL2Cfoam.set_OMP(shell_parallelization)
    verbosity_flux && myid() == 1 && log("sl2cfoam-next initialized on each worker with C parallelization set to $(shell_parallelization)")

end