# preliminary check
function check_on_preliminary_parameters(data_folder_path::String, sl2cfoam_next_data_folder::String, Dl_min::Int64, Dl_max::Int64, number_of_workers::Int, number_of_threads::Int)

    if (number_of_workers * number_of_threads > length(Sys.cpu_info()))
        printstyled("WARNING: you are using more resources than available cores on this system. Performances will be affected\n\n"; bold=true, color=:red)
    end
    sleep(1)

    # data folder check
    if (isdir(data_folder_path) == false)
        error("data folder path does not exists")
    end

    # sl2cfoam_next data folder check
    if (isdir(sl2cfoam_next_data_folder) == false)
        error("data folder path of sl2cfoam_next does not exists")
    end

    # shell parameter check
    if (typeof(Dl_min) != Int64 || Dl_min < 0)
        error("Please assign number of minimum shell as positive integer")
    end
    if (typeof(Dl_max) != Int64 || Dl_max < 0)
        error("Please assign number of maximum shell as positive integer")
    end
    if (Dl_max < Dl_min)
        error("Maximum number of shell must be larger than the minimum one")
    end

end


# check each input configuration
function check_configuration!(user_conf)

    j = user_conf[1]

    # set spin to float (ending in .0) if integer
    if (typeof(j) == Int64)
        user_conf[1] = convert(Float64, user_conf[1])
        j = convert(Float64, j)
    end

    # spin check
    if (j < 0)
        error("Please assign positive spin in user_conf $(user_conf)\n")
    end

    #TODO: add further checks

end