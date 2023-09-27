using LinearAlgebra
using Zygote
include("transfer_map_zygote.jl")
include("linear_optics.jl")

# function FODO_track_result(k1, k2)
#     r = [0.01,0.001,0.0,0.0,0.0,0.0]
#     QF = FQuad("QF", 1, k1,0.0) 
#     D1 = Drift("D1", 1)
#     QD = DFQuad("QD", 1, k2, 0.0)
#     D2 = Drift("D2",1)
#     FODO = [QF,D1,QD,D2] 

#     r1 = track(FODO, [0.01,0.001,0.0,0.0,0.0,0.0])
#     y = sqrt(r1[1]^2 + r1[2]^2 + r1[3]^2 + r1[4]^2 + r1[5]^2 + r1[6]^2)
#     return y
# end

# k1 = -3.0
# k2 = 2.0
# # FODO_track_result(k1,k2)

# # Compute the gradient
# gradient_g = Zygote.gradient(FODO_track_result, k1, k2) #error Mutating arrays is not supported
# print(gradient_g)


function f(k1, k2, L1, L2)
    QF = Quad("QF", 1, k1, 0.0)
    D1 = Drift("D1", L1)
    QD = Quad("QD", 1, k2, 0.0)
    D2 = Drift("D2",L1)
    FODO = [QF,D1,QD,D2]
    alpha0 = [0, 0]
    beta0 = [3, 5]
    alphax, alphy, betax, betay, mux, muy = twissline(FODO, 0, alpha0, beta0)
    return betax[1] + betay[1]
end

using Zygote
k1 = -3
k2 = 2
L1 = 1
L2 = 2
# gradient_g = Zygote.gradient(f, k1, k2)
gradient_g = gradient(params -> f(params[1], params[2], L1, L2), [k1, k2])
print("gradient at k1 = $k1, k2 = $k2: $gradient_g.")

# function f(k1, k2)
#     QF = Quad("QF", 1, k1, 0.0)
#     D1 = Drift("D1", 1)
#     QD = Quad("QD", 1, k2, 0.0)
#     D2 = Drift("D2",1)
#     FODO = [QF,D1,QD,D2]
#     rin = [0.01 0.005;
#             0.001 0.002;
#             0 0.01;
#             0 0.001;
#             0 0;
#             0 0]
#     refpts = [1, 2, 3, 5]
#     rout = linepass1(FODO, rin, refpts)
#     return rout
# end

# k1 = -3
# k2 = 2
# rout = f(k1,k2)
# print(rout)

