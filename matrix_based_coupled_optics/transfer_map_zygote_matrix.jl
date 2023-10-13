# module elements
include("beam.jl")
using LinearAlgebra
using Zygote
# export Drift, FQuad, DFQuad, Bend, TransferMap, track

##############################################

abstract type AbstractElement end
struct Drift <: AbstractElement
    name::String
    len::Float64
end

struct Quad <: AbstractElement
    name::String
    len::Float64
    k1::Float64 
    k1s::Float64
end

struct ThinQuad <: AbstractElement
    name::String
    len::Float64
    k1l::Float64 
    k1sl::Float64
    function ThinQuad(name::String, k1l::Float64, k1sl::Float64)
        return new(name, 0.0, k1l, k1sl)
    end
end

struct SBend <: AbstractElement
    name::String
    len::Float64
    angle::Float64
    k1::Float64
    e1::Float64
    e2::Float64
    f_int1::Float64
    f_int2::Float64
end

struct RBend <: AbstractElement
    name::String
    len::Float64
    angle::Float64
    k1::Float64
    e1::Float64
    e2::Float64
    f_int1::Float64
    f_int2::Float64
end

struct DipEdge <: AbstractElement
    name::String
    len::Float64
    h::Float64
    e1::Float64
    f_int::Float64
    function DipEdge(name::String, h::Float64, e1::Float64, f_int::Float64)
        return new(name, 0.0, h, e1, f_int)
    end
end

struct DipBody <: AbstractElement
    name::String
    len::Float64
    angle::Float64
    k1::Float64
end

struct ThinCrabCavity <: AbstractElement
    name::String
    len::Float64
    strength::Tuple{Float64,Float64}
    frequency::Float64
    phase::Float64
    kcc::Float64
end
function ThinCrabCavity(;name::String="",strength::Tuple{Float64,Float64},frequency::Float64,phase::Float64=0.0)
	kcc::Float64=2pi*Frequency/Float64(299792458)
	ThinCrabCavity(name,0,strength,frequency,phase,kcc)
end

struct Marker <: AbstractElement
    name::String
    len::Float64
    function Marker(name::String)
        return new(name, 0.0)
    end
end

struct Solenoid <: AbstractElement
    name::String
    len::Float64
    ks::Float64
end

struct LorentzBoost <: AbstractElement
    name::String
    len::Float64
    thetac::Float64
    thetas::Float64
    function LorentzBoost(name::String,thetac::Float64,thetas::Float64)
        return new(name, 0.0, thetac, thetas)
    end
end

struct RevLorentzBoost <: AbstractElement
    name::String
    len::Float64
    thetac::Float64
    thetas::Float64
    function RevLorentzBoost(name::String,thetac::Float64,thetas::Float64)
        return new(name, 0.0, thetac, thetas)
    end
end
##############################################
# Drift
function transferMatrix(mag::Drift,beam::AbstractBeam=_beam[])
    M = [1.0 mag.len 0 0 0 0;
        0 1.0 0 0 0 0;
        0 0 1.0 mag.len 0 0;
        0 0 0 1.0 0 0;
        0 0 0 0 mag.len/(beam.betagamma0^2) 0;
        0 0 0 0 0 1.0]
    return M	
end

# Bend
function _edgeMatrix(theta::Float64,h0::Float64)
	t = tan(theta)*h0
    M = [1.0 0 0 0 0 0;
        t 1.0 0 0 0 0;
        0 0 1.0 0 0 0;
        0 0 -t 1.0 0 0;
        0 0 0 0 1.0 0;
        0 0 0 0 0 1.0]
	return M
end

function _sectorBodyMatrix(L::Float64, h0::Float64, dx::Float64, 
                            Kx::Float64, Ky::Float64, rx::Float64, ry::Float64, beta0::Float64, betagamma0::Float64)
    M11, M12, M13, M14, M15, M16 = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    M21, M22, M23, M24, M25, M26 = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    M31, M32, M33, M34, M35, M36 = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    M41, M42, M43, M44, M45, M46 = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    M51, M52, M53, M54, M55, M56 = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    M61, M62, M63, M64, M65, M66 = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0

    if rx > Float64(0)
        D0 = dx / Kx
        tx = sqrt(Kx)
        th = tx * L
        C, S = cos(th), sin(th)

        M11, M22 = C, C
        M12 = S / tx
        M21 = -tx * S
        M16 = D0 * (1.0 - C)
        M26 = D0 * S * tx
        M56 = -D0 * (L * h0 - S * h0 / tx)

    elseif rx < Float64(0)
        D0 = dx / Kx
        tx = sqrt(-Kx)
        th = tx * L
        C, S = cosh(th), sinh(th)

        M11, M22 = C, C
        M12 = S / tx
        M21 = tx * S
        M16 = D0 * (1 - C)
        M26 = -D0 * S * tx
        M56 = -D0 * (L * h0 - S * h0 / tx)

    else
        M11, M22 = 1.0, 1.0
        M12 = L
        M21 = 0.0
        M16 = L * L * dx / 2
        M26 = dx * L
        M56 = -dx * h0 * L^3 / 6
    end

    if ry > Float64(0)
        ty = sqrt(Ky)
        th = ty * L
        C, S = cos(th), sin(th)

        M33, M44 = C, C
        M34 = S / ty
        M43 = -ty * S

    elseif ry < Float64(0)
        ty = sqrt(-Ky)
        th = ty * L
        C, S = cosh(th), sinh(th)

        M33, M44 = C, C
        M34 = S / ty
        M43 = ty * S

    else
        M33, M44 = 1.0, 1.0
        M34 = L
    end

    M16 /= beta0
    M26 /= beta0

    M56 /= beta0 * beta0
    M56 += L / (betagamma0^2)

    M15 = M12 * M11 - M11 * M22
    M25 = M12 * M21 - M11 * M26

    M55, M66 = 1.0, 1.0

    return [
        M11 M12 M13 M14 M15 M16;
        M21 M22 M23 M24 M25 M26;
        M31 M32 M33 M34 M35 M36;
        M41 M42 M43 M44 M45 M46;
        M51 M52 M53 M54 M55 M56;
        M61 M62 M63 M64 M65 M66
    ]
end

function transferMatrix(mag::SBend,beam::AbstractBeam=_beam[])
	mag.len == Float64(0) && (return Float64(1)*Matrix(I,(6,6)))
	h0 = mag.angle/mag.len
	Kx = h0*h0+mag.k1
	Ky = -mag.k1
	rx = Kx/(h0*h0);ry=Ky/(h0*h0)
	dx = h0

	M = _sectorBodyMatrix(mag.len,h0,dx,Kx,Ky,rx,ry,beam.beta0,beam.betagamma0)

	if mag.e1 != Float64(0)
		M = M*_edgeMatrix(mag.e1,h0)
	end

	if mag.e2 != Float64(0)
		M = _edgeMatrix(mag.e2,h0)*M
	end

	return M
end

function transferMatrix(mag::RBend,beam::AbstractBeam=_beam[])
	mag.len == Float64(0) && (return Float64(1)*Matrix(I,(6,6)))
	h0 = Float64(2)*sin(mag.angle/Float64(2))/mag.len
	Kx = h0*h0+mag.k1
	Ky = -mag.k1
	rx = Kx/(h0*h0);ry=Ky/(h0*h0)
	dx = h0

	M = _sectorBodyMatrix(mag.angle/h0,h0,dx,Kx,Ky,rx,ry,beam.beta0,beam.betagamma0)

	return _edgeMatrix(mag.e2+mag.angle/Float64(2),h0)*M*_edgeMatrix(mag.e1+mag.angle/Float64(2),h0)
end

function transferMatrix(mag::DipEdge,beam::AbstractBeam=_beam[])
	return _edgeMatrix(mag.e1,mag.h)
end

function transferMatrix(mag::DipBody,beam::AbstractBeam=_beam[])
	mag.len == Float64(0) && (return Float64(1)*Matrix(I,(6,6)))
	h0 = mag.angle/mag.len
	Kx = h0*h0+mag.k1
	Ky = -mag.k1
	rx = Kx/(h0*h0);ry=Ky/(h0*h0)
	dx = h0
	return _sectorBodyMatrix(mag.len,h0,dx,Kx,Ky,rx,ry,beam.beta0,beam.betagamma0)
end

# Quadrupole
function transferMatrix(mag::Quad, beam::AbstractBeam = _beam[])
    b = mag.k1
    a = mag.k1s
    t1 = sqrt(b * b + a * a)

    if t1 == Float64(0)
        return [
            1.0   mag.len   0.0   0.0   0.0   0.0;
            0.0   1.0   0.0   0.0   0.0   0.0;
            0.0   0.0   1.0   mag.len   0.0   0.0;
            0.0   0.0   0.0   1.0   0.0   0.0;
            0.0   0.0   0.0   0.0   1.0   mag.len/(beam.betagamma0^2);
            0.0   0.0   0.0   0.0   0.0   1.0
        ]
    end

    C, S = b + t1, -a
    t2 = sqrt(C * C + S * S)
    if t2 != Float64(0)
        C /= t2
        S /= t2
    else
        C, S = -a, b - t1
        t2 = sqrt(C * C + S * S)
        C /= t2
        S /= t2
    end
    t1 = sqrt(t1)
    th = t1 * mag.len

    A = [
        C   0.0  S   0.0;
        0.0 C   0.0 S;
        -S  0.0 C   0.0;
        0.0 -S  0.0 C
    ]

    invA = [
        C   0.0  -S   0.0;
        0.0 C   0.0 -S;
        S  0.0 C   0.0;
        0.0 S  0.0 C
    ]
    
    M44 = [cos(th) sin(th)/t1 0.0 0.0;
        -sin(th)*t1 cos(th) 0.0 0.0;
        0.0 0.0 cosh(th) sinh(th)/t1;
        0.0 0.0 sinh(th)*t1 cosh(th)]
    
    M44 = invA * M44 * A
    # Build the 6*6 map matrix using M44 as a submatrix 
    M = [M44 zeros(Float64, 4, 2); 0.0 0.0 0.0 0.0 1.0 mag.len/(beam.betagamma0^2); 0.0 0.0 0.0 0.0 0.0 1.0]
    
    return M
end

function transferMatrix(mag::ThinQuad, beam::AbstractBeam = _beam[])
    M = [
        1.0      0.0      0.0      0.0      0.0      0.0;
        -mag.k1l 1.0      mag.k1sl 0.0      0.0      0.0;
        0.0      0.0      1.0      0.0      0.0      0.0;
        mag.k1sl 0.0      mag.k1l  1.0      0.0      0.0;
        0.0      0.0      0.0      0.0      1.0      0.0;
        0.0      0.0      0.0      0.0      0.0      1.0
    ]
    return M
end

# Solenoid
function transferMatrix(mag::Solenoid, beam::AbstractBeam = _beam[])
    K = mag.ks/Float64(2)
    th = K * mag.len
    C, S = cos(th), sin(th)
    CC = C * C
    SC = S * C
    SS = S * S
    ret = [
        CC       SC/K      SC      SS/K      0.0      0.0;
        -SC*K      CC       -K*SS      SC      0.0      0.0;
        -SC      -SS/K      CC       SC/K      0.0      0.0;
        K*SS      -SC      -SC*K      CC       0.0      0.0;
        0.0      0.0      0.0      0.0      mag.len/(beam.betagamma0^2)      0.0;
        0.0      0.0      0.0      0.0      0.0      1.0
    ]
    return ret
end

# Crab cavity
function transferMatrix(mag::ThinCrabCavity,beam::AbstractBeam=_beam[])
    M=[
        1.0      0.0      0.0      0.0      0.0      0.0;
        0.0      1.0      0.0      0.0      -mag.strength[1]      0.0;
        0.0      0.0      1.0      0.0      0.0      0.0;
        0.0      0.0      0.0      1.0      -mag.strength[2]      0.0;
        0.0      0.0      0.0      0.0      1.0      0.0;
        -mag.strength[1]      0.0      -mag.strength[2]      0.0      0.0      1.0
    ]
    return M
end

# Boost
function transferMatrix(mag::LorentzBoost,beam::AbstractBeam=_beam[])
    M = [
        1.0      0.0      0.0      0.0      Float64(1)/cos(mag.thetac)*sin(mag.thetac)      0.0;
        0.0      Float64(1)/cos(mag.thetac)      0.0      0.0      0.0      0.0;
        0.0      0.0      1.0      0.0      0.0      0.0;
        0.0      0.0      0.0      Float64(1)/cos(mag.thetac)      0.0      0.0;
        0.0      0.0      0.0      0.0      Float64(1)/cos(mag.thetac)      0.0;
        0.0      -Float64(1)/cos(mag.thetac)*sin(mag.thetac)      0.0      0.0      0.0      1.0
    ]
    if mag.thetas==Float64(0)
        return M
    else
        R = [
        cos(mag.thetas)      0.0      sin(mag.thetas)      0.0      0.0      0.0;
        0.0      cos(mag.thetas)      0.0      sin(mag.thetas)      0.0      0.0;
        -sin(mag.thetas)      0.0      cos(mag.thetas)      0.0      0.0      0.0;
        0.0      -sin(mag.thetas)      0.0      cos(mag.thetas)      0.0      0.0;
        0.0      0.0      0.0      0.0      1.0      0.0;
        0.0      0.0      0.0      0.0      0.0      1.0
    ]
        return M*R 
    end 
end

function transferMatrix(mag::RevLorentzBoost,beam::AbstractBeam=_beam[])
    M = [
        1.0      0.0      0.0      0.0      -sin(mag.thetac)      0.0;
        0.0      cos(mag.thetac)      0.0      0.0      0.0      0.0;
        0.0      0.0      1.0      0.0      0.0      0.0;
        0.0      0.0      0.0      cos(mag.thetac)      0.0      0.0;
        0.0      0.0      0.0      0.0      cos(mag.thetac)      0.0;
        0.0      sin(mag.thetac)      0.0      0.0      0.0      1.0
    ]
    if mag.thetas==Float64(0)
        return M
    else
        R = [
        cos(mag.thetas)      0.0      -sin(mag.thetas)      0.0      0.0      0.0;
        0.0      cos(mag.thetas)      0.0      -sin(mag.thetas)      0.0      0.0;
        sin(mag.thetas)      0.0      cos(mag.thetas)      0.0      0.0      0.0;
        0.0      sin(mag.thetas)      0.0      cos(mag.thetas)      0.0      0.0;
        0.0      0.0      0.0      0.0      1.0      0.0;
        0.0      0.0      0.0      0.0      0.0      1.0
    ]
        return R*M 
    end 
end

# Marker
transferMatrix(mag::Marker,beam::AbstractBeam=_beam[])=Float64(1)*Matrix(I,(6,6))




# struct Particle
#     E::Number
#     r::Vector{Float64}
# end

# function track(cell, r::Vector{Float64})
#     r1 = copy(r)
#     for ele in cell
#         r1, map = TransferMap(ele, r1)
#     end
#     return r1
# end

# function linepass(cell, rin)
#     np = length(rin[1,:])
#     r_buff = Zygote.Buffer(rin)
#     for i in 1:np
#         rin_m = [rin[1,i], rin[2,i]* (1 + rin[5]), rin[3,i], rin[4,i]* (1 + rin[5]), rin[5,i], rin[6,i]]
#         rout = track(cell, rin_m)
#         r_buff[:,i] = rout
#     end
#     return copy(r_buff)
# end

# # end
# function track_through_refpt(cell, r::Vector{Float64}, refpts::Vector{Int})
#     r1 = copy(r)
#     results = Zygote.Buffer(r, length(r), length(refpts))

#     ref_index = 1
#     current_refpt = refpts[ref_index]

#     # If the first reference point is 1, store the initial state
#     if current_refpt == 1
#         results[:, ref_index] = r1
#         ref_index += 1
#         if ref_index <= length(refpts)
#             current_refpt = refpts[ref_index]
#         else
#             return copy(results)
#         end
#     end

#     for (i, ele) in enumerate(cell)
#         r1, _ = TransferMap(ele, r1)

#         # Check if current element is a reference point
#         if i == current_refpt - 1
#             results[:, ref_index] = r1

#             ref_index += 1
#             if ref_index <= length(refpts)
#                 current_refpt = refpts[ref_index]
#             else
#                 break
#             end
#         end
#     end

#     return copy(results)
# end

# function linepass1(cell, rin, refpts=nothing)
#     np = size(rin, 2)
#     r_buff = Zygote.Buffer(rin, size(rin, 1), np * (refpts === nothing ? 1 : length(refpts)))


#     for i in 1:np
#         rin_m = [rin[1,i], rin[2,i]* (1 + rin[5,i]), rin[3,i], rin[4,i]* (1 + rin[5,i]), rin[5,i], rin[6,i]]

#         # If refpts is specified, store the state at each refpt.
#         if refpts !== nothing
#             rout = track_through_refpt(cell, rin_m, refpts)
#             r_buff[:, (i-1)*size(rout, 2) .+ (1:size(rout, 2))] = rout
#         else
#             rout = track(cell, rin_m)
#             r_buff[:,i] = rout
#         end
#     end

#     return copy(r_buff)
# end



