include("transfer_map_zygote_derong.jl")
include("beam.jl")
using StaticArrays
using StructArrays
abstract type AbstractTwiss end

struct EdwardsTengTwiss <: AbstractTwiss
	betax::Float64
	betay::Float64
	alphax::Float64
	alphay::Float64
	gammax::Float64
	gammay::Float64
	dx::Float64
	dpx::Float64
	dy::Float64
	dpy::Float64
	sinmux::Float64
	cosmux::Float64
	sinmuy::Float64
	cosmuy::Float64
	# alpha::SVector{2,Float64}
	# gamma::SVector{2,Float64}
	# eta::SVector{4,Float64}
	# sc_mu::SVector{4,Float64} # sin(mu1) cos(mu1) sin(mu2) cos(mu2)
	R::Matrix{Float64}
	mode::Int
end

EdwardsTengTwiss(;betx::Float64,bety::Float64,
			   alfx::Float64=Float64(0),alfy::Float64=Float64(0),
			   dx::Float64=Float64(0),dy::Float64=Float64(0),
			   dpx::Float64=Float64(0),dpy::Float64=Float64(0),
			   mux::Float64=Float64(0),muy::Float64=Float64(0),
			   R11::Float64=Float64(0),R12::Float64=Float64(0),
			   R21::Float64=Float64(0),R22::Float64=Float64(0),
			   mode::Int=Int(1))=EdwardsTengTwiss(betx,bety,alfx,alfy,(Float64(1)+alfx^2)/betx,(Float64(1)+alfy^2)/bety,
														  dx,dpx,dy,dpy,sin(mux),cos(mux),sin(muy),cos(muy),[R11 R12;R21 R22],mode)
_symplectic_conjugate_2by2(M) = [M[2, 2] -M[1, 2]; -M[2, 1] M[1, 1]]
_matrixTransform_2by2(M)=begin
    m11,m21,m12,m22=M
    [m11*m11 -2m11*m12 m12*m12
    -m11*m21 1.0+2m12*m21 -m12*m22
    m21*m21 -2m21*m22 m22*m22]
end

function twissPropagate(tin::EdwardsTengTwiss,M::Matrix{Float64})
	A=@view M[1:2,1:2]
	B=@view M[1:2,3:4]
	C=@view M[3:4,1:2]
	D=@view M[3:4,3:4]

	R1=tin.R
	_R1=_symplectic_conjugate_2by2(R1)
	if tin.mode == Int(1)
		X=A-B*R1
		begin t=det(X)
			if t>Float64(0.1)
				R=(D*R1-C)*_symplectic_conjugate_2by2(X)
				R/=t
				X/=sqrt(t)
				Y=D+C*_R1
				Y/=sqrt(det(Y))
				mode=Int(1)
			else
				X=C-D*R1
				X/=sqrt(det(X))
				Y=B+A*_R1
				t=det(Y)
				R=-(D+C*_R1)*_symplectic_conjugate_2by2(Y)
				R/=t
				Y/=sqrt(t)
				mode=Int(2)
			end
		end
	elseif tin.mode == Int(2) 
		X=B+A*_R1
		begin t=det(X)
			if t>Float64(0.1)
				R=-(D+C*_R1)*_symplectic_conjugate_2by2(X)
				R/=t
				X/=sqrt(t)
				Y=C-D*R1
				Y/=sqrt(det(Y))
				mode=Int(1)
			else
				X=D+C*_R1
				X/=sqrt(det(X))
				Y=A-B*R1
				t=det(Y)
				R=(D*R1-C)*_symplectic_conjugate_2by2(Y)
				R/=t
				Y/=sqrt(t)
				mode=Int(2)
			end
		end
	else
		#throw(AssertionError("Mode should be integer 1 or 2."))
		println(stderr,"Invalid mode.")
		return EdwardsTengTwiss(;betx=Float64(1),bety=Float64(1),mode=Int(0))
	end

	Nx=_matrixTransform_2by2(X)
	Ny=_matrixTransform_2by2(Y)
	v1=Nx*[tin.betax;tin.alphax;tin.gammax]
	v2=Ny*[tin.betay;tin.alphay;tin.gammay]
	eta=(@view M[1:4,1:4])*[tin.dx,tin.dpx,tin.dy,tin.dpy]+(@view M[1:4,6])
	sin_dmux=X[1,2]/sqrt(v1[1]*tin.betax)
	cos_dmux=X[1,1]*sqrt(tin.betax/v1[1])-tin.alphax*sin_dmux
	sin_dmuy=Y[1,2]/sqrt(v2[1]*tin.betay)
	cos_dmuy=Y[1,1]*sqrt(tin.betay/v2[1])-tin.alphay*sin_dmuy

	smux0,cmux0,smuy0,cmuy0=tin.sinmux,tin.cosmux,tin.sinmuy,tin.cosmuy
	smux=sin_dmux*cmux0+cos_dmux*smux0
	cmux=cos_dmux*cmux0-sin_dmux*smux0
	smuy=sin_dmuy*cmuy0+cos_dmuy*smuy0
	cmuy=cos_dmuy*cmuy0-sin_dmuy*smuy0

	return EdwardsTengTwiss(v1[1],v2[1],v1[2],v2[2],v1[3],v2[3],eta[1],eta[2],eta[3],eta[4],smux,cmux,smuy,cmuy,R,mode)
end

function twissPropagate(tin::EdwardsTengTwiss,seq::Vector{AbstractElement},beam::AbstractBeam=_beam[])
	# this is not good for AD
	ret=Vector{typeof(tin)}(undef,1+length(seq))
	ss=zeros(Float64,length(ret))
	names=Vector{String}(undef,length(ret))
	ret[1]=tin
	names[1]="Start"
	for (index,mag) in enumerate(seq)
		M=transferMatrix(mag,beam)
		ret[index+1]=twissPropagate(ret[index],M)
		ss[index+1]=mag.len + ss[index]
		names[index+1]=mag.name
	end
	return ss,names,StructArray(ret)
end

function twissPropagate(tin::EdwardsTengTwiss,seq::Vector{AbstractElement},beam::AbstractBeam=_beam[], endindex::Int=1)
    ret = tin
    ss = 0.0
	for (index,mag) in enumerate(seq)
		M=transferMatrix(mag,beam)
		ret=twissPropagate(ret,M)
		ss=mag.len + ss
		names=mag.name
		if index == endindex
			break
		end
	end
	return ss,names,ret
end
function periodicEdwardsTengTwiss(M::Matrix{Float64})
	A=@view M[1:2,1:2]
	B=@view M[1:2,3:4]
	C=@view M[3:4,1:2]
	D=@view M[3:4,3:4]
	invalid_ret=EdwardsTengTwiss(;betx=Float64(1),bety=Float64(1),mode=Int(0))

	Bbar_and_C=_symplectic_conjugate_2by2(B)+C
	t1=0.5*(tr(A)-tr(D))
	Δ=t1*t1+det(Bbar_and_C)
	Δ<Float64(0) && (println(stderr,"Failed to decouple periodic transfer matrix. The linear matrix is unstable.");return invalid_ret)

	_sign= t1>Float64(0) ? Float64(-1) : Float64(1)

	t2=abs(t1)+sqrt(Δ)
	if t2==Float64(0)
		R=Float64[0 0;0 0]
	else
		R=Bbar_and_C*(_sign/t2)
	end

	X=A-B*R
	Y=D+C*_symplectic_conjugate_2by2(R)

	# It should be equal to 1
	(det(X)<Float64(0.9) || det(Y)<Float64(0.9))  && (println(stderr,"Failed to decouple the periodic transfer matrix with mode 1.");return invalid_ret)

	cmux=Float64(0.5)*(X[1,1]+X[2,2])
	cmuy=Float64(0.5)*(Y[1,1]+Y[2,2])
	(Float64(-1)<cmux<Float64(1) && Float64(-1)<cmuy<Float64(1)) || (println(stderr,"Failed to get beta functions. The linear matrix is unstable.");return invalid_ret)

	smux=sqrt(Float64(1)-cmux*cmux)*sign(X[1,2])
	smuy=sqrt(Float64(1)-cmuy*cmuy)*sign(Y[1,2])
	betx=X[1,2]/smux
	gamx=-X[2,1]/smux
	bety=Y[1,2]/smuy
	gamy=-Y[2,1]/smuy

	alfx=Float64(0.5)*(X[1,1]-X[2,2])/smux
	alfy=Float64(0.5)*(Y[1,1]-Y[2,2])/smuy

	eta=inv(Matrix{Float64}(I,(4,4))-(@view M[1:4,1:4]))*(@view M[1:4,6])
	return EdwardsTengTwiss(betx,bety,alfx,alfy,gamx,gamy,eta[1],eta[2],eta[3],eta[4],smux,cmux,smuy,cmuy,R,Int(1))
end


function normalMatrix(tin::EdwardsTengTwiss)
	(tin.mode==Int(1) || tin.mode==Int(2)) || begin
		println(stderr,"Warning: return identity matrix for unknown mode $(tin.mode) as the normal matrix (transformation matrix from normal space to physical space).")
		return Float64(1)*Matrix{Float64}(I,6,6)
	end
	D=Float64[1 0 0 0 0 tin.eta[1]
			   0 1 0 0 0 tin.eta[2]
			   0 0 1 0 0 tin.eta[3]
			   0 0 0 1 0 tin.eta[4]
			   -tin.eta[2] tin.eta[1] -tin.eta[4] tin.eta[3] 1 0
			   0 0 0 0 0 1]
	sbx,sby=sqrt.(tin.beta)
	B=Float64[sbx 0 0 0 0 0
			   -tin.alphax/sbx 1/sbx 0 0 0 0
			   0 0 sby 0 0 0
			   0 0 -tin.alphay/sby 1/sby 0 0
			   0 0 0 0 1 0
			   0 0 0 0 0 1]
	λ=Float64(1)/sqrt(abs(Float64(1)+det(tin.R)))
	R=λ*tin.R
	_R=_symplectic_conjugate_2by2(R)
	O=Float64[0 0;0 0]
	U=Float64[λ 0;0 λ]
	if tin.mode==Int(1)
		V=[U _R O;-R U O;O O I]
	else
		V=[_R U O;U -R O;O O I]
	end
	return D*V*B
end

