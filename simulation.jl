
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
