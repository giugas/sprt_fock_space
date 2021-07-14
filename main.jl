cd(dirname(@__FILE__))
using QuantumOptics
using LinearAlgebra
using Plots
using JLD

#using DifferentialEquations (⊗ is defined in both QuantumOptics and DifferentialEquations)
include("functions.jl")
include("function_innovation.jl")

# Parameters
ω_mech = 100.
ω_opt=10^4
Δ = ω_opt-ω_mech
# Constants
g = 1.
η = 1.
k = 0.5
γ  = 0.1
n= 10
rates =[ k 0 0; 0 γ*(n+1) 0; 0 0 γ*n ]

# Basis
b_cav = FockBasis(20,0)
b_mech = FockBasis(20,0)


# Operators Cavity
a = destroy(b_cav)⊗one(b_mech)
at = create(b_cav) ⊗ one(b_mech)

# Operators Oscillator
b = one(b_cav) ⊗ destroy(b_mech)
bt = one(b_cav) ⊗ create(b_mech)

# Hamilton operator
H_cav(Δ,η) = Δ*at*a + η*(a - at)
H_mech(ω_mech) = -ω_mech*bt*b
H_int(g) = -g*(bt+b)*at*a

H_tot(Δ,ω_mech,g,η) = H_cav(Δ,η) + H_mech(ω_mech) + H_int(g)


H=H_tot(Δ,ω_mech,g,η)
# Diffusive term operators
J     = [a, b,bt]
Jdagger =dagger.(J)

#measurmtent Operator
θ=  π/2
C=-sqrt(k)im*a
#C=-sqrt(k)*a*exp(-1.0im *θ)
Cdagger = dagger(C)
dt=0.00001


# Initial State of the system
ψ = (fockstate(b_cav,20)⊗fockstate(b_mech,10))#+fockstate(b_cav,9)⊗fockstate(b_mech,20))/sqrt(2.)
ψ = coherentstate(b_cav, 0.0001)⊗coherentstate(b_mech, 0.01+0.0001*im)
ρ = ψ⊗dagger(ψ)




N = 5000 #number of steps
dt = 0.000002 #step pace

Y=time_trace(N,dt,ρ,C,J)
#save("/tmp/myfile.jld", "time_trace", Y)

#Y = load("/tmp/myfile.jld","time_trace")
ts= 0:dt:(N-1)*dt
plot(real.(Y[2,1:3000]))
plot!(real.(Y[1,1:3000]))


function simulation(
    n::Int,
    Y::Vector,
    dt::Float64,
    h1,
    h0,
    )

    H1=h1[1]
    rates1=h1[2]
    η1 = h1[3]
    k1 = h1[4]
    ρ1 = h1[5]

    H0=h0[1]
    rates0=h0[2]
    η0 = h0[3]
    k0 = h0[4]
    ρ0 = h0[5]

    tmp1 = copy(ρ)
    tmp0 = copy(ρ)
    dρ = copy(ρ)
    Y=real(Y)

    L=Array{Float64}(undef,4,n + 1)
    L[:,1] = [0,0,0,0]
    for i = 2:n
        x1=real(tr(χ(ρ1,C,Cdagger,tmp1)))
        L[1,i] = L[1,i-1]+loglike(η1,x1,Y[i-1])
        L[2,i]=x1
        ρ1+=dmaster(ρ1,H1,rates1,J,Jdagger,dρ,tmp1)*dt+
        dmaster_stochastic(ρ1,C,Cdagger,tmp1)*innovation(Y[i-1],η1,x1,dt)

        x0=real(tr(χ(ρ0,C,Cdagger,tmp0)))
        L[3,i] = L[3,i-1]+loglike(η0,x0,Y[i-1])
        L[4,i]= x0
        ρ0+=dmaster(ρ0,H0,rates0,J,Jdagger,dρ,tmp0)*dt+
        dmaster_stochastic(ρ0,C,Cdagger,tmp1)*innovation(Y[i-1],η0,x0,dt)

        #L[i]=1-l0
        print("Processing step $i\u001b[1000D")
    end
    return L
end




ω = ω_mech
Δ = ω_opt-ω_mech
ρ1 = ψ⊗dagger(ψ)
ρ0 = ψ⊗dagger(ψ)
H1 = H_tot(Δ,ω,g,η)
H0 = H_tot(Δ,ω,g,η)
H  = H_tot(Δ,ω,g,η)
rates1=rates
rates0=rates
C1=C
C1dagger=Cdagger
C0=C
C0dagger=Cdagger
J1=J
J0=J
J1dagger=Jdagger
J0dagger=Jdagger
η1=η
η0=η
k1=k
k0=k

h1=(H1,rates1, 0.5*η1 , k1 ,ρ1 )
h0=(H0,rates0, η0 , k0 ,ρ0 )

no = 3000
Y1 = real(Y[2,1:no])
L = simulation(no,Y1,dt,h1,h0)
plot(L[1,1:3000]-L[3,1:3000])
plot(L[4,1:2000])
plot(Y1[:])
