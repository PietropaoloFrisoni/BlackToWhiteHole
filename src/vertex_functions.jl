# distributes the computation of many vertex tensors over workers on a single machine
function VertexDistributeSingleMachine(spins_configurations, Dl::Integer, store=(true, true))

    @sync @distributed for spins in spins_configurations
        vertex_compute(collect(spins), Dl; result=VertexResult((false, store[1], store[2])))
    end

end

# computes a coherent state coefficient for the BW spinfoam boundary
# each direction is +1 if TARGET, -1 if SOURCE
# angles are received as [[theta1, phi1], [theta2, phi2] ... ]
function CoherentStateVector!(coherentstate::Vector{ComplexF64}, spins::Vector{HalfInt8}, angles::Vector{Vector{Float64}}, directions::Vector{Int64})

    j1 = spins[1]
    j2 = spins[2]
    j3 = spins[3]
    j4 = spins[4]

    theta1 = angles[1][1]
    phi1 = angles[1][2]

    theta2 = angles[2][1]
    phi2 = angles[2][2]

    theta3 = angles[3][1]
    phi3 = angles[3][2]

    theta4 = angles[4][1]
    phi4 = angles[4][2]

    sgn1 = half(2 * directions[1])
    sgn2 = half(2 * directions[2])
    sgn3 = half(2 * directions[3])
    sgn4 = half(2 * directions[4])

    if (sgn1 == 1)
        theta1 -= pi / 2
    end
    if (sgn2 == 1)
        theta2 -= pi / 2
    end
    if (sgn3 == 1)
        theta3 -= pi / 2
    end
    if (sgn4 == 1)
        theta4 -= pi / 2
    end

    range_tuple = intertwiner_range(j1, j2, j3, j4)

    for intertwiner_index = 1:range_tuple[2]

        intertwiner = half(2 * from_index_to_intertwiner(range_tuple, intertwiner_index))

        for m1::HalfInt8 = -j1:j1
            for m2::HalfInt8 = -j2:j2

                if (m1 + m2 > intertwiner || m1 + m2 < -intertwiner)
                    continue
                end

                for m3::HalfInt8 = -j3:j3

                    m4 = -m1 - m2 - m3

                    if (m4 > j4 || m4 < -j4)
                        continue
                    end

                    W4j = Wigner4jm(j1, j2, j3, j4, m1, m2, m3, m4, intertwiner)

                    WignerD1_with_sign = wignerDjmn(j1, m1, j1 * sgn1, phi1, theta1, 0.0)
                    WignerD2_with_sign = wignerDjmn(j2, m2, j2 * sgn2, phi2, theta2, 0.0)
                    WignerD3_with_sign = wignerDjmn(j3, m3, j3 * sgn3, phi3, theta3, 0.0)
                    WignerD4_with_sign = wignerDjmn(j4, m4, j4 * sgn4, phi4, theta4, 0.0)

                    if (sgn1 == -1)
                        WignerD1_with_sign *= (-1)^(j1 + m1)
                    end
                    if (sgn2 == -1)
                        WignerD2_with_sign *= (-1)^(j2 + m2)
                    end
                    if (sgn3 == -1)
                        WignerD3_with_sign *= (-1)^(j3 + m3)
                    end
                    if (sgn4 == -1)
                        WignerD4_with_sign *= (-1)^(j4 + m4)
                    end

                    @inbounds coherentstate[intertwiner_index] += W4j * WignerD1_with_sign * WignerD2_with_sign * WignerD3_with_sign * WignerD4_with_sign

                end
            end
        end

        @inbounds coherentstate[intertwiner_index] *= sqrt(Dimension(intertwiner))

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

        counter += 1

        i2_range = intertwiners_range[index][2]
        i3_range = intertwiners_range[index][3]
        i4_range = intertwiners_range[index][4]
        i5_range = intertwiners_range[index][5]

        vertex = vertex_load([spins_configurations[index][1], spins_configurations[index][2], spins_configurations[index][3],
                spins_configurations[index][4], spins_configurations[index][5], spins_configurations[index][6],
                spins_configurations[index][7], spins_configurations[index][8], spins_configurations[index][9],
                spins_configurations[index][10]], Dl)

        coherentstate1 = Vector{ComplexF64}(undef, i2_range)
        coherentstate2 = Vector{ComplexF64}(undef, i3_range)
        coherentstate3 = Vector{ComplexF64}(undef, i4_range)
        coherentstate4 = Vector{ComplexF64}(undef, i5_range)

        coherentstate1[:] .= 0.0 + 0.0 * im
        coherentstate2[:] .= 0.0 + 0.0 * im
        coherentstate3[:] .= 0.0 + 0.0 * im
        coherentstate4[:] .= 0.0 + 0.0 * im

        CoherentStateVector!(coherentstate1,
            [spins_configurations[index][5], spins_configurations[index][6], spins_configurations[index][7], spins_configurations[index][1]],
            [[theta, phi_1], [theta, phi_2], [theta, phi_3], [0.0, 0.0]], [1, 1, -1, 1])

        CoherentStateVector!(coherentstate2,
            [spins_configurations[index][8], spins_configurations[index][9], spins_configurations[index][2], spins_configurations[index][5]],
            [[theta, phi_3], [theta, phi_2], [0.0, 0.0], [theta, phi_1]], [-1, 1, -1, -1])

        CoherentStateVector!(coherentstate3,
            [spins_configurations[index][10], spins_configurations[index][3], spins_configurations[index][6], spins_configurations[index][8]],
            [[theta, phi_2], [0.0, 0.0], [theta, phi_3], [theta, phi_1]], [1, 1, -1, 1])

        CoherentStateVector!(coherentstate4,
            [spins_configurations[index][4], spins_configurations[index][7], spins_configurations[index][9], spins_configurations[index][10]],
            [[0.0, 0.0], [theta, phi_2], [theta, phi_1], [theta, phi_3]], [-1, 1, -1, -1])

        for i1 in 1:i1_range
            s = 0.0 + 0.0 * im
            @turbo warn_check_args = false for i2 in 1:i2_range, i3 in 1:i3_range, i4 in 1:i4_range, i5 in 1:i5_range
                s += vertex.a[i5, i4, i3, i2, i1] * coherentstate1[i2] * coherentstate2[i3] * coherentstate3[i4] * coherentstate4[i5]
            end

            coherent_matrix_up[i1, counter] = s

        end

    end

end

# contracts the vertex tensor with the lower boundary coherent state
# DEPRECATED
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
            [[theta, 0, theta, theta] [phi_2, 0, phi_3, phi_1]])

        c4 = coherentstate_compute([spins_configurations[index][4]
                spins_configurations[index][7]
                spins_configurations[index][9]
                spins_configurations[index][10]],
            [[0, theta, theta, theta] [0, phi_2, phi_1, phi_3]])

        counter += 1

        extra_rotation_factor = exp(-im * (0 * spins_configurations[index][5] + (phi_2 - phi_3) * spins_configurations[index][6] + (phi_2 - phi_3) * spins_configurations[index][7]
                                           + (phi_1 - phi_3) * spins_configurations[index][8] + (phi_2 - phi_1) * spins_configurations[index][9] + (phi_2 - phi_3) * spins_configurations[index][10]))

        for i1 in 1:i1_range
            s = 0.0 + 0.0 * im
            @turbo warn_check_args = false for i2 in 1:i2_range, i3 in 1:i3_range, i4 in 1:i4_range, i5 in 1:i5_range
                s += vertex.a[i5, i4, i3, i2, i1] * c1.a[i2] * c2.a[i3] * c3.a[i4] * c4.a[i5]
            end

            coherent_matrix_down[i1, counter] = s * extra_rotation_factor

        end

    end

end

# contracts the upper and lower coherent matrices
function SpinfoamFinalContraction!(contracted_spinfoam, coherent_matrix_up, coherent_matrix_down, i1_range, total_radial_spins_combinations)

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

end