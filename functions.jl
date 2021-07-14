function dmaster(rho::T, H::AbstractOperator{B,B},
                    rates::Matrix, J::Vector, Jdagger::Vector,
                    drho::T, tmp::T) where {B<:Basis,T<:Operator{B,B}}
    QuantumOpticsBase.mul!(drho,H,rho,-eltype(rho)(im),zero(eltype(rho)))
    QuantumOpticsBase.mul!(drho,rho,H,eltype(rho)(im),one(eltype(rho)))
    for j=1:length(J), i=1:length(J)
        QuantumOpticsBase.mul!(tmp,J[i],rho,eltype(ρ)(rates[i,j]),zero(eltype(ρ)))
        QuantumOpticsBase.mul!(drho,tmp,Jdagger[j],true,true)
        QuantumOpticsBase.mul!(drho,Jdagger[j],tmp,eltype(rho)(-0.5),one(eltype(rho)))
        QuantumOpticsBase.mul!(tmp,rho,Jdagger[j],eltype(rho)(rates[i,j]),zero(eltype(rho)))
        QuantumOpticsBase.mul!(drho,tmp,J[i],eltype(rho)(-0.5),one(eltype(rho)))
    end
    return drho

end

function dmaster(rho::T, H::AbstractOperator{B,B},
                    rates::Float64, J::AbstractOperator{B,B}, Jdagger::AbstractOperator{B,B},
                    drho::T, tmp::T) where {B<:Basis,T<:Operator{B,B}}
    QuantumOpticsBase.mul!(drho,H,rho,-eltype(rho)(im),zero(eltype(rho)))
    QuantumOpticsBase.mul!(drho,rho,H,eltype(rho)(im),one(eltype(rho)))
    QuantumOpticsBase.mul!(tmp,J,rho,eltype(ρ)(rates),zero(eltype(ρ)))
    QuantumOpticsBase.mul!(drho,tmp,Jdagger,true,true)
    QuantumOpticsBase.mul!(drho,Jdagger,tmp,eltype(rho)(-0.5),one(eltype(rho)))
    QuantumOpticsBase.mul!(tmp,rho,Jdagger,eltype(rho)(rates),zero(eltype(rho)))
    QuantumOpticsBase.mul!(drho,tmp,J,eltype(rho)(-0.5),one(eltype(rho)))
    return drho
end


#dmastper(rho)= Cρ+ρCdagger-Tr((C+Cdagger)ρ)
function dmaster_stochastic(rho::T,
            C::AbstractOperator{B,B}, Cdagger::AbstractOperator{B,B}, drho::T) where {B<:Basis,T<:Operator{B,B}}
        QuantumOpticsBase.mul!(drho,C,rho,false,false)
        QuantumOpticsBase.mul!(drho,rho,Cdagger,true,true)
        drho-= tr(drho)*rho
    return drho
end



#... X = tr(χ(ρ,C,Cdagger,tmp))
# dw=innovation(dI,η,X,dt)
# drho=dmaster(rho,H,rates,J,Jdagger,drho,tmp)*dt+dmaster_stochastic(rho,C,Cdagger,drho)*dw
