# distributes the computation of many vertex tensors over workers on a single machine
function VertexDistributeSingleMachine(spins_configurations, Dl::Integer, store=(true, true))

    @sync @distributed for spins in spins_configurations
        vertex_compute(collect(spins), Dl; result=VertexResult((false, store[1], store[2])))
    end

end

# contracts the vertex tensor with the upper boundary coherent state
function SpinfoamContractUp!(coherent_matrix_up, lower_bound, upper_bound, i1_range, spins_configurations, Dl)

    theta = acos(-sqrt(2 / 3))
    phi_1 = 0
    phi_2 = 2 * pi / 3
    phi_3 = -2 * pi / 3

    counter = 0

    @inbounds for index = lower_bound:upper_bound 

        i2_range = intertwiners_range[index][2]
        i3_range = intertwiners_range[index][3]
        i4_range = intertwiners_range[index][4]
        i5_range = intertwiners_range[index][5]

        vertex = vertex_load([spins_configurations[index][1], spins_configurations[index][2], spins_configurations[index][3],
                spins_configurations[index][4], spins_configurations[index][5], spins_configurations[index][6],
                spins_configurations[index][7], spins_configurations[index][8], spins_configurations[index][9],
                spins_configurations[index][10]], Dl)

        c1 = coherentstate_compute([spins_configurations[index][5]
                spins_configurations[index][6]
                spins_configurations[index][7]
                spins_configurations[index][1]],
            [[theta, theta, theta, 0] [phi_1, phi_2, phi_3, 0]])

        c2 = coherentstate_compute([spins_configurations[index][8]
                spins_configurations[index][9]
                spins_configurations[index][2]
                spins_configurations[index][5]],
            [[theta, theta, 0, theta] [phi_3, phi_2, 0, phi_1]])

        c3 = coherentstate_compute([spins_configurations[index][10]
                spins_configurations[index][3]
                spins_configurations[index][6]
                spins_configurations[index][8]],
            [[theta, 0, theta, theta] [phi_1, 0, phi_2, phi_3]])

        c4 = coherentstate_compute([spins_configurations[index][4]
                spins_configurations[index][7]
                spins_configurations[index][9]
                spins_configurations[index][10]],
            [[0, theta, theta, theta] [0, phi_3, phi_2, phi_1]])

        counter += 1

        for i1 in 1:i1_range
            s = 0.0 + 0.0 * im
            @turbo warn_check_args = false for i2 in 1:i2_range, i3 in 1:i3_range, i4 in 1:i4_range, i5 in 1:i5_range
                s += vertex.a[i5, i4, i3, i2, i1] * c1.a[i2] * c2.a[i3] * c3.a[i4] * c4.a[i5]
            end

            coherent_matrix_up[i1, counter] = s

        end

    end

end

# contracts the vertex tensor with the lower boundary coherent state
function SpinfoamContractDown!(coherent_matrix_down, lower_bound, upper_bound, i1_range, spins_configurations, Dl)

    theta = acos(-sqrt(2 / 3))
    phi_1 = 0
    phi_2 = 2 * pi / 3
    phi_3 = -2 * pi / 3

    counter = 0

    @inbounds for index = lower_bound:upper_bound

        i2_range = intertwiners_range[index][2]
        i3_range = intertwiners_range[index][3]
        i4_range = intertwiners_range[index][4]
        i5_range = intertwiners_range[index][5]

        vertex = vertex_load([spins_configurations[index][1], spins_configurations[index][2], spins_configurations[index][3],
                spins_configurations[index][4], spins_configurations[index][5], spins_configurations[index][6],
                spins_configurations[index][7], spins_configurations[index][8], spins_configurations[index][9],
                spins_configurations[index][10]], Dl)

        c1 = coherentstate_compute([spins_configurations[index][5]
                spins_configurations[index][6]
                spins_configurations[index][7]
                spins_configurations[index][1]],
            [[theta, theta, theta, 0] [phi_1, phi_2, phi_3, 0]])

        c2 = coherentstate_compute([spins_configurations[index][8]
                spins_configurations[index][9]
                spins_configurations[index][2]
                spins_configurations[index][5]],
            [[theta, theta, 0, theta] [phi_3, phi_2, 0, phi_1]])

        c3 = coherentstate_compute([spins_configurations[index][10]
                spins_configurations[index][3]
                spins_configurations[index][6]
                spins_configurations[index][8]],
            [[theta, 0, theta, theta] [phi_1, 0, phi_2, phi_3]])

        c4 = coherentstate_compute([spins_configurations[index][4]
                spins_configurations[index][7]
                spins_configurations[index][9]
                spins_configurations[index][10]],
            [[0, theta, theta, theta] [0, phi_3, phi_2, phi_1]])

        counter += 1

        for i1 in 1:i1_range
            s = 0.0 + 0.0 * im
            @turbo warn_check_args = false for i2 in 1:i2_range, i3 in 1:i3_range, i4 in 1:i4_range, i5 in 1:i5_range 
                s += vertex.a[i5, i4, i3, i2, i1] * c1.a[i2] * c2.a[i3] * c3.a[i4] * c4.a[i5]
            end

            coherent_matrix_down[i1, counter] = s

        end

    end

end

# contracts the upper and lower coherent matrices
function SpinfoamFinalContraction!(contracted_spinfoam, coherent_matrix_up, coherent_matrix_down, i1_range, total_radial_spins_combinations, j1, j2, j3, j4, Dl, immirzi)

    counter = 0

    @inbounds for index1 = 1:total_radial_spins_combinations   

        @inbounds for index2 = 1:total_radial_spins_combinations

            s = 0.0 + 0.0 * im

            @turbo warn_check_args = false for i1 in 1:i1_range
                s += coherent_matrix_up[i1, index1] * coherent_matrix_down[i1, index2]
            end

            counter += 1
            contracted_spinfoam[counter] = s

        end
        
    end

    @save "$(conf.spinfoam_folder)/j1_$(j1)_j2_$(j2)_j3_$(j3)_j4_$(j4)/immirzi_$(immirzi)/contracted_spinfoam_Dl_$(Dl).jld2" contracted_spinfoam

end