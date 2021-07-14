



# dρ = Cρ+ρCdagger
function χ(rho::T,C::AbstractOperator{B,B},Cdagger::AbstractOperator{B,B},tmp::T) where {B<:Basis,T<:Operator{B,B}}
    QuantumOpticsBase.mul!(tmp,C,rho)
    QuantumOpticsBase.mul!(tmp,rho,Cdagger,true,true)
    return tmp
end

#dW innovation term!!
#X= Tr(χ)

function innovation(dI::Float64,η::Float64,X::Float64,dt::Float64)
    return  dI-η*X*dt
end

#(the dI= ηXφdt+dW)
function measured_signal(X::Float64,η::Float64,noise::Float64,dt::Float64)
    return η*X*dt+noise
end

#loglike for a single process.
#loglike = √η *X*dI-0.5*η*X^2*dI²  if dI²∝dt and where X=tr(C+Cdagger rho) .
function loglike(η::Float64,x::Float64,dI::Float64)
    tmp= sqrt(η)*x*dI
    dl =  (1.0-0.5*tmp)*tmp
    return dl
end

function loglike(η::Float64,x::Float64,dI::Float64,dt)
    tmp= sqrt(ηk)*dI
    dl =  x*tmp -0.5*η*x^2*dt
    return dl
end



function time_trace(n::Int , dt::Float64, ρ::T,
        C::AbstractOperator{B,B},J::Vector) where {B<:Basis,T<:Operator{B,B}}
    Y=Array{ComplexF64}(undef,2,n )
    tmp = copy(ρ)
    dρ = copy(ρ)
    Cdagger=dagger(C)
    for i = 1:n
        dw= randn()*sqrt(dt)
        Y[1,i] = tr(χ(ρ,C,Cdagger,tmp))
        Y[2,i] = measured_signal(Y[1,i].re,η,dw,dt)
        ρ += dmaster(ρ,H,rates,J,Jdagger,dρ,tmp)*dt+sqrt(η)*dmaster_stochastic(ρ,C,Cdagger,dρ)*dw

        print("Processing step $i\u001b[1000D")
    end
    return Y
end

function time_trace(n::Int , dt::Float64, ρ::T, C::AbstractOperator{B,B},J::Vector) where {B<:Basis,T<:Operator{B,B}}
    Y=Array{ComplexF64}(undef,2,n )
    tmp = copy(ρ)
    dρ = copy(ρ)
    Cdagger=dagger(C)
    for i = 1:n
        dw= randn()*sqrt(dt)
        ρ+= dmaster(ρ,H,rates,J,Jdagger,dρ,tmp)*dt+sqrt(η)*dmaster_stochastic(ρ,C,Cdagger,dρ)*dw
        Y[1,i] = tr(χ(ρ,C,Cdagger,tmp))
        Y[2,i] =  measured_signal(real(Y[1,i]), η , dw , dt )
        print("Processing step $i\u001b[1000D")
    end
    return Y
end
