# preliminary check
function CheckPreliminaryParameters(data_folder_path::String, sl2cfoam_next_data_folder::String, Dl_min::Int64, Dl_max::Int64, number_of_workers::Int, number_of_threads::Int)

    # resources check
    (number_of_workers * number_of_threads > length(Sys.cpu_info())) && warn("more procs than available cores on this system. Performances will be affected")

    # data folder check
    !isdir(data_folder_path) && error("data folder path does not exists")

    # sl2cfoam_next data folder check
    !isdir(sl2cfoam_next_data_folder) && error("data folder path of sl2cfoam_next does not exists")

    # truncation parameter checks
    Dl_min < 0 && error("Assign positive or null minimum truncation parameter")
    Dl_max < 0 && error("Assign positive or null maximum truncation parameter")
    Dl_max < Dl_min && error("Maximum value of the truncation parameter cannot be smaller than minimum")

    # remaining checks
    immirzi < 0 && error("Assign positive Immirzi parameter")
    T_sampling_parameter <= 0 && error("Assign positive T sampling parameter")
    size(angular_spins, 1) == 0 && error("Assign positive number of spins configurations")

end


# check each input configuration
function CheckConfiguration!(user_conf)

    j = user_conf[1]

    # set spin to float (ending in .0) if integer
    if (typeof(j) == Int64)
        user_conf[1] = convert(Float64, user_conf[1])
        j = convert(Float64, j)
    end

    # spin check
    j < 0 && error("Please assign positive spin in user_conf $(user_conf)\n")

end

# check if there's NaN in (multidimensional) array
function sum_check(x)
    s = sum(x)
    isnan(s) || !isfinite(s)
end