mutable struct Configuration

    j0::Half{Int8}
    j0_float::Float64
    K0::Half{Int8}
    jpm::Half{Int8}
    jpm_float::Float64
    Kpm::Half{Int8}

    contracted_spinfoam_found::Bool

    base_folder::String
    spinfoam_folder::String

end

function init_config(user_conf, data_folder_path::String, contracted_spinfoam_found::Bool=false)

    j0_float = user_conf[1]
    jpm_float = round(Int64, 2 * (j0_float / (sqrt(6)))) / 2
    j0 = half(2 * j0_float)
    jpm = half(2 * jpm_float)

    K0 = half(2)
    Kpm = half(2)

    # ensure that there's at least 1 possible value (triangular inequalities)
    if (j0_float == 7.5)
        jpm_float = 3.5
        jpm = half(2 * jpm_float)
    end
    if (j0_float == 8)
        jpm_float = 3
        jpm = half(2 * jpm_float)
    end
    if (j0_float == 9)
        jpm_float = 4
        jpm = half(2 * jpm_float)
    end
    if (j0_float == 9.5)
        jpm_float = 3.5
        jpm = half(2 * jpm_float)
    end
    if (j0_float == 11)
        jpm_float = 4
        jpm = half(2 * jpm_float)
    end
    if (j0_float == 12.5)
        jpm_float = 5.5
        jpm = half(2 * jpm_float)
    end

    intermediate_path_data_folder = "/data/spinfoam_data/j0=$(j0_float)_jpm=$(jpm_float)"

    # where spins_conf, spins_map and intertwiner_range are stored
    base_folder = "$(data_folder_path)/$(intermediate_path_data_folder)"

    # where ??? is stored
    spinfoam_folder = "$(base_folder)/spinfoam"

    conf = Configuration(j0, j0_float, K0, jpm, jpm_float, Kpm,
        contracted_spinfoam_found,
        base_folder, spinfoam_folder)

    return conf

end

function init_sl2cfoam_next(Immirzi_parameter, sl2cfoam_next_data_folder, number_of_threads, verbosity_flux)

    isMPI = @ccall SL2Cfoam.clib.sl2cfoam_is_MPI()::Bool
    isMPI && error("MPI version not allowed")

    conf_sl2cfoam_next = SL2Cfoam.Config(VerbosityOff, HighAccuracy, 100, 0)
    SL2Cfoam.cinit(sl2cfoam_next_data_folder, Immirzi_parameter, conf_sl2cfoam_next)

    shell_parallelization = false

    # enable C library automatic parallelization
    if (number_of_threads > 1)
        shell_parallelization = true
    end

    SL2Cfoam.set_OMP(shell_parallelization)

    if (verbosity_flux && Bool(myid()))
        println("sl2cfoam-next initialized on each worker with C parallelization set to $(shell_parallelization)")
    end

end