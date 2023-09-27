# module linear_optics
# This optics calculation code is the same as AT, except the settings of optimization iterations
using LinearAlgebra
include("transfer_map_zygote.jl")

function findm44(LATTICE, dp=0.0)
    NE = length(LATTICE)
    XYStep = 3e-8  # Default step size for numerical differentiation
    scaling = XYStep * ones(4)
    
    # Call findorbit4 function to find initial orbit
    # Define findorbit4 function according to your needs
    _, orbitin = findorbit4(LATTICE, dp)

    # Add end-of-lattice
    refs = vcat(1:NE, NE+1)

    # Build a diagonal matrix of initial conditions
    D4 = [0.5 * Diagonal(scaling);zeros(2,4)]
    
    # Add to the orbit_in. First 8 columns for derivative
    # 9-th column is for closed orbit
    RIN = [orbitin .+ D4 orbitin .- D4 orbitin]

    # Call linepass function to propagate particles
    # Define linepass function according to your needs
    ROUT = linepass(LATTICE, RIN)
    
    TMAT3 = reshape(ROUT[1:4,:], (4, 9, :))
    M44 = (TMAT3[:,1:4,end] - TMAT3[:,5:8,end]) ./ scaling

    return M44
end

function findspos(line, refpts)
    L = [0; cumsum(map(el -> el.len, line))]
    return L[refpts]
end

function findorbit4(ring, dp)
    # Default orbit initial condition is [0; 0; 0; 0; 0; 0]
    orbitin = xorbit_dp(ring, dp)
    orb4 = orbitin[1:4, 1]
    return orb4, orbitin
end

function xorbit_dp(ring, dp=0.0)
    XYStep = 3e-8   # Step size for numerical differentiation
    dps = 1e-12    # Convergence threshold
    max_iterations = 20  # Max. iterations
    Ri = reshape([0; 0; 0; 0; dp; 0], 6, 1)  # Default initial guess for the orbit

    scaling = XYStep .* [1, 1, 1, 1]
    D = [Diagonal(scaling) zeros(4,1); zeros(2,5)]

    # change = Inf

    for itercount in 1:max_iterations
        RMATi = Ri * ones(1,5) + D
        RMATf = linepass(ring, RMATi)
        Rf = RMATf[:, end]
        # Compute the transverse part of the Jacobian
        J4 = (RMATf[1:4,1:4] .- RMATf[1:4,5]) ./ scaling
        Ri_next = Ri + [ (I - J4) \ (Rf[1:4] - Ri[1:4]); 0; 0]
        # change = norm(Ri_next - Ri)
        Ri = Ri_next
        # if change < dps
        #     break
        # end
    end

    return Ri
end

function twissline(LINE, DP, alpha, beta)
    NE = length(LINE)

    ax = alpha[1]
    ay = alpha[2]
    bx = beta[1]
    by = beta[2]

    # M44, MS, orb = findm44(LINE, DP, NE + 1, R0)
    MS = findm44(LINE, DP)
    BX = ((MS[1,1,:] .* bx .- MS[1,2,:] .* ax) .^ 2 .+ MS[1,2,:].^2) ./ bx
    BY = ((MS[3,3,:] .* by .- MS[3,4,:] .* ay) .^ 2 .+ MS[3,4,:].^2) ./ by
    AX = - ((MS[1,1,:] .* bx .- MS[1,2,:] .* ax) .* (MS[2,1,:] .* bx .- MS[2,2,:] .* ax) .+ MS[1,2,:] .* MS[2,2,:]) ./ bx
    AY = - ((MS[3,3,:] .* by .- MS[3,4,:] .* ay) .* (MS[4,3,:] .* by .- MS[4,4,:] .* ay) .+ MS[3,4,:] .* MS[4,4,:]) ./ by
    MX = BetatronPhaseUnwrap(atan.(MS[1,2,:] ./ (MS[1,1,:] .* bx .- MS[1,2,:] .* ax)))
    MY = BetatronPhaseUnwrap(atan.(MS[3,4,:] ./ (MS[3,3,:] .* by .- MS[3,4,:] .* ay)))
    # TD = Dict(
    #     "ElemIndex" => NE + 1,
    #     "SPos" => findspos(LINE, NE + 1),
    #    # "ClosedOrbit" => orb,
    #     "M44" => MS,
    #     "beta" => hcat(BX,BY),
    #     "alpha" => hcat(AX,AY),
    #     "mu" => hcat(MX,MY)
    # )
    return AX, AY, BX, BY, MX, MY
end

function BetatronPhaseUnwrap(P)
    DP = diff(P)
    JUMPS = vcat(0, DP .< -1e-3)
    return P + cumsum(JUMPS) * pi
end

# end

