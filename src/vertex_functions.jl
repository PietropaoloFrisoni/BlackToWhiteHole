# distributes the computation of many vertex tensors over workers on a single machine
function vertex_distributed_single_machine(spins_configurations, Dl::Integer, store=(true, true))

    @sync @distributed for spins in spins_configurations

        vertex_compute(collect(spins), Dl; result=VertexResult((false, store[1], store[2])))
        
    end

end