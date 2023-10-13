abstract type AbstractBeam end
struct Beam <: AbstractBeam 
	E0::Float64
	m0::Float64
	gamma0::Float64
	beta0::Float64
	betagamma0::Float64
	Brho::Float64
end

Beam(;E0::Float64, m0::Float64)=begin
	gamma=E0/m0
	betagamma=sqrt(gamma*gamma-1)
	beta=betagamma/gamma
	Brho=sqrt(E0^2-m0^2)*Float64(1e9)/Float64(299792458)
	Beam(E0,m0,gamma,beta,betagamma,Brho)
end

const _beam=(Beam(E0=Float64(275.0),m0=Float64(0.93827208816)) |> Ref)

setDefaultBeam(;E0::Float64, m0::Float64)=begin
	global _beam[]=Beam(;E0=E0,m0=m0)
	return
end

showDefaultBeam()=show(_beam[])

macro DefaultBeam() _beam[] end
