module Comparison

const KS_MNPA = include("KS_MNPA.jl").KS_MNPA
const BFF_MNPA = include("BFF_MNPA.jl").BFF_MNPA
const KS_Pauli = include("KS_pauli.jl").KS_Pauli
const KS_Mel = include("KS_mel.jl").KS_Mel

function comparison()
    q = 0.3

    println("Running Matrix BFF method")
    res = @time BFF_MNPA.vN_entropy_BB84(q)
    @show res
    GC.gc()

    println("Running Matrix KS method")
    res = @time KS_MNPA.vN_entropy_BB84(q)
    @show res
    GC.gc()

    println("Running Pauli KS method")
    res = @time KS_Pauli.vN_entropy(q)
    @show res
    GC.gc()

    println("Running Matrix KS method")
    res = @time KS_Mel.vN_entropy(q)
    @show res
    GC.gc()
end

end # module Comparison
