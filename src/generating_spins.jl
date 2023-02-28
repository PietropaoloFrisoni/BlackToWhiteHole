# generate a list of all spins to compute and make folders
function generating_spins(j0, K0, jpm, Kpm, base_folder, spinfoam_folder, immirzi, verbosity_flux)

    spins_configurations = NTuple{10,HalfInt8}[]
    intertwiners_range = NTuple{5,Int8}[]
    onehalf = half(1)
  
    #number_of_angular_spins_combinations = file_count(vertex_folder)
    #spins_map = Vector{Int64}(undef, number_of_angular_spins_combinations)
  
    spins_map = ElasticArray{Int64}(undef, 0)
  
    counter_angular_spins_combinations = 0
  
    for j12::HalfInt8 = (j0-K0):onehalf:(j0+K0),
      j13::HalfInt8 = (j0-K0):onehalf:(j0+K0),
      j14::HalfInt8 = (j0-K0):onehalf:(j0+K0),
      j15::HalfInt8 = (j0-K0):onehalf:(j0+K0)
  
      # skip if any intertwiner range empty
      r1, r1_range = intertwiner_range(j12, j13, j14, j15)
      isempty(r1) && continue
  
      mkpath(spinfoam_folder * "/j1_$(twice(j12)/2)_j2_$(twice(j13)/2)_j3_$(twice(j14)/2)_j4_$(twice(j15)/2)/immirzi_$(immirzi)")
  
      counter_angular_spins_combinations += 1
  
      number_of_radial_spins_combinations_for_this_angular_spin_combination = 0
  
      for j23::HalfInt8 = (jpm-Kpm):onehalf:(jpm+Kpm),
        j24::HalfInt8 = (jpm-Kpm):onehalf:(jpm+Kpm),
        j25::HalfInt8 = (jpm-Kpm):onehalf:(jpm+Kpm),
        j34::HalfInt8 = (jpm-Kpm):onehalf:(jpm+Kpm),
        j35::HalfInt8 = (jpm-Kpm):onehalf:(jpm+Kpm),
        j45::HalfInt8 = (jpm-Kpm):onehalf:(jpm+Kpm)
  
        # skip if any intertwiner range empty
        r2, r2_range = intertwiner_range(j12, j25, j24, j23)
        r3, r3_range = intertwiner_range(j23, j13, j34, j35)
        r4, r4_range = intertwiner_range(j34, j24, j14, j45)
        r5, r5_range = intertwiner_range(j45, j35, j25, j15)
  
        isempty(r2) && continue
        isempty(r3) && continue
        isempty(r4) && continue
        isempty(r5) && continue
  
        number_of_radial_spins_combinations_for_this_angular_spin_combination += 1
  
        # must be computed
        push!(spins_configurations, (j12, j13, j14, j15, j23, j24, j25, j34, j35, j45))
        push!(intertwiners_range, (r1_range, r2_range, r3_range, r4_range, r5_range))
  
      end
  
      resize!(spins_map, counter_angular_spins_combinations)
      spins_map[counter_angular_spins_combinations] = number_of_radial_spins_combinations_for_this_angular_spin_combination
  
    end
  
    if (verbosity_flux == true)
  
      println("Number of angular spins combinations: $(counter_angular_spins_combinations)\n",
        "Total number of spins combinations: $(size(spins_configurations)[1])\n")
  
    end
  
    @save "$(base_folder)/spins_configurations.jld2" spins_configurations
    @save "$(base_folder)/intertwiners_range.jld2" intertwiners_range
    @save "$(base_folder)/spins_map.jld2" spins_map
  
  end