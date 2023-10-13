include("EdwardsTengTwiss.jl")
include("Beam.jl")
include("transfer_map_zygote_matrix.jl")
include("sequence.jl")
include("RipkenTwiss.jl")
using Zygote

function f(k1, k2, L1, L2)
    D1 = Drift("D1", L1)
    D2 = Drift("D2", L2)
    QD = Quad("QD", 1, k2, 0)
    QF = Quad("QF", 1, k1, 0)
    twissin = EdwardsTengTwiss(betx=1.0, bety=1.0, alfx=0.0, alfy=0.0, dx=0.0, dy=0.0, dpx=0.0, 
                                dpy=0.0, mux=0.0, muy=0.0, R11=0.0, R12=0.0, R21=0.0, R22=0.0, mode=1)
    seq = [QF, D1, QD, D2]
    beam = Beam(E0=275.0, m0=0.93827208816)
    ss, elementname, twissout = twissPropagate(twissin, seq, beam, 4)
    betax = twissout.betax
    betay = twissout.betay
    return betax + betay
end
function f1(k1, k2, L1, L2)
    D1 = Drift("D1", L1)
    D2 = Drift("D2", L2)
    QD = Quad("QD", 1, k2, 0)
    QF = Quad("QF", 1, k1, 0)
    twissin = RipkenTwiss(1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0,
                            0.0, 1.0, 0.0, 1.0, 0.0)
    seq = [QF, D1, QD, D2]
    beam = Beam(E0=275.0, m0=0.93827208816)
    ss, elementname, twissout = twissPropagate(twissin, seq, beam, 2)
    betax = twissout.betaxx
    betay = twissout.betayy
    return betax + betay
end

k1 = -3
k2 = 2
L1 = 0.5
L2 = 0.5
init = f(k1, k2, L1, L2)
print(init)
# @time gradient_g = gradient(params -> f(params[1], params[2], L1, L2), [k1, k2])
@time gradient_g = gradient(params -> f(params[1], params[2], L1, L2), [k1, k2])
print("gradient at k1 = $k1, k2 = $k2: $gradient_g.")
perturb = f(k1 + 0.0001, k2, L1, L2)
delta = (perturb - init) / 0.0001
print(delta)