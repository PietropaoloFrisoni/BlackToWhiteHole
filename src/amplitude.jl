function WeightFactor!(weight_factor, alpha, j0, jpm, m, T_range, Immirzi, spins_configurations,
    lower_bound, upper_bound, j1, j2, j3, j4)

    zita_plus = -32 * sqrt(6) / 9
    zita_minus = 32 * sqrt(6) / 9

    pre_fact_angul_1 = exp(-(j0^alpha) / (2))
    pre_fact_rad = exp(-(jpm^alpha) / (2))

    pre_fact_plus = exp(im * Immirzi * zita_plus)
    pre_fact_minus = exp(im * Immirzi * zita_minus)

    @inbounds for T_index in eachindex(T_range)

        zita_0 = T_range[T_index] / (2 * m)
        pre_fact_angul_2 = exp(im * Immirzi * zita_0)
        angular_factor = (2 * j1 + 1) * (2 * j2 + 1) * (2 * j3 + 1) * (2 * j4 + 1) *
                         pre_fact_angul_2^(j1) * pre_fact_angul_2^(j2) * pre_fact_angul_2^(j3) * pre_fact_angul_2^(j4) *
                         pre_fact_angul_1^((j1 - j0)^2) * pre_fact_angul_1^((j2 - j0)^2) *
                         pre_fact_angul_1^((j3 - j0)^2) * pre_fact_angul_1^((j4 - j0)^2)

        weight_factor[:, T_index] .= 0.0 + 0.0 * im

        counter = 0

        @inbounds for index_minus = lower_bound:upper_bound

            j12_minus = spins_configurations[index_minus][5]
            j13_minus = spins_configurations[index_minus][6]
            j14_minus = spins_configurations[index_minus][7]
            j23_minus = spins_configurations[index_minus][8]
            j24_minus = spins_configurations[index_minus][9]
            j34_minus = spins_configurations[index_minus][10]

            minus_radial_factor = (2 * j12_minus + 1) * (2 * j13_minus + 1) * (2 * j14_minus + 1) *
                                  (2 * j23_minus + 1) * (2 * j24_minus + 1) * (2 * j34_minus + 1) *
                                  pre_fact_rad^((j12_minus - jpm)^2) * pre_fact_rad^((j13_minus - jpm)^2) * pre_fact_rad^((j14_minus - jpm)^2) *
                                  pre_fact_rad^((j23_minus - jpm)^2) * pre_fact_rad^((j24_minus - jpm)^2) * pre_fact_rad^((j34_minus - jpm)^2) *
                                  pre_fact_minus^(j12_minus) * pre_fact_minus^(j13_minus) * pre_fact_minus^(j14_minus) *
                                  pre_fact_minus^(j23_minus) * pre_fact_minus^(j24_minus) * pre_fact_minus^(j34_minus)

            @inbounds for index_plus = lower_bound:upper_bound

                counter += 1

                j34_plus = spins_configurations[index_plus][5]
                j24_plus = spins_configurations[index_plus][6]
                j14_plus = spins_configurations[index_plus][7]
                j23_plus = spins_configurations[index_plus][8]
                j13_plus = spins_configurations[index_plus][9]
                j12_plus = spins_configurations[index_plus][10]

                plus_radial_factor = (2 * j34_plus + 1) * (2 * j24_plus + 1) * (2 * j14_plus + 1) *
                                     (2 * j23_plus + 1) * (2 * j13_plus + 1) * (2 * j12_plus + 1) *
                                     pre_fact_rad^((j34_plus - jpm)^2) * pre_fact_rad^((j24_plus - jpm)^2) * pre_fact_rad^((j14_plus - jpm)^2) *
                                     pre_fact_rad^((j23_plus - jpm)^2) * pre_fact_rad^((j13_plus - jpm)^2) * pre_fact_rad^((j12_plus - jpm)^2) *
                                     pre_fact_plus^(j34_plus) * pre_fact_plus^(j24_plus) * pre_fact_plus^(j14_plus) *
                                     pre_fact_plus^(j23_plus) * pre_fact_plus^(j13_plus) * pre_fact_plus^(j12_plus)

                radial_factor = minus_radial_factor * plus_radial_factor

                weight_factor[counter, T_index] += radial_factor * angular_factor

            end

        end

    end

end

function AmplitudeIntegration!(amplitude_abs_sq_integrated, amplitude_abs_sq, T_range)

    Delta_x = T_range[2] - T_range[1]
    amp = 0.0
    #amp_T = 0.0

    for Dl_index in eachindex(amplitude_abs_sq_integrated)

        for T_index in eachindex(T_range)

            if (T_index == 1)
                amp += final_data_vector_abs_sq[T_index]
                amp_T += final_data_vector_abs_sq_T[T_index]
            end

            amp += 2 * final_data_vector_abs_sq[T_index]
            amp_T += 2 * final_data_vector_abs_sq_T[T_index]

            if (T_index == conf.total_T_points)
                amp += final_data_vector_abs_sq[T_index]
                amp_T += final_data_vector_abs_sq_T[T_index]
            end

        end

        Integrated_amplitude = (Delta_x / 2) * amp
        Integrated_amplitude_T = (Delta_x / 2) * amp_T

        amplitude_abs_sq_integrated[:]

    end

end