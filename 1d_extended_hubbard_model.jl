## DMRG calculation of the one-dimensional extended Hubbard model
## Used to compare to ED result
using ITensors

let 
    N = 20
    Npart = 20
    t1 = 1.0             # Nearest-neighbor hopping amplitude
    t2 = 0.0             # Next-nearest-neighbor hopping amplitude
    U = 1.0              # On-site Coulomb interaction
    V = 0.5              # Nearest-neighbor interaction 



    sites = siteinds("Electron", N; conserve_qns = true)

    # Set up the 1d extended Hubbard Hamiltonian
    ampo = OpSum()

    # Nearest-neighbor terms
    # Using open boundary condition
    for bond in 1 : (N - 1)
        ampo += -t1, "Cdagup", bond, "Cup", bond + 1
        ampo += -t1, "Cdagup", bond + 1, "Cup", bond
        ampo += -t1, "Cdagdn", bond, "Cdn", bond + 1
        ampo += -t1, "Cdagdn", bond + 1, "Cdn", bond
        ampo += V1, "Ntot", bond, "Ntot", bond + 1
    end

    # Next-nearest-neighbor terms
    # Using open boundary condition
    for bond in 1 : (N - 2)
        ampo += -t2, "Cdagup", bond, "Cup", bond + 2
        ampo += -t2, "Cdagup", bond + 2, "Cup", bond
        ampo += -t2, "Cdagdn", bond, "Cdn", bond + 2
        ampo += -t2, "Cdagdn", bond + 2, "Cdn", bond
    end

    for ind in 1 : N
        ampo += U, "Nupdn", ind
    end

    H = MPO(ampo, sites)

    nsweeps = 10
    maxdim = [50, 100, 200, 400, 800, 1000]
    cutoff = [1E-12]

    state = ["Emp" for n in 1 : N]
    tmp = Npart
    for ind in N - 1:-1:1
        if tmp > ind
            println("Doubly occupying site $ind")
            state[ind] = "UpDn"
            tmp -= 2
        elseif tmp > 0
            println("Singly occupying site $ind")
            state[ind] = (isodd(ind) ? "Up" : "Dn")
            tmp -= 1
        end
    end 

    # Initialize wavefunction using random MPS
    ψ₀ = randomMPS(sites, state, 10)

    @show flux(ψ₀)

    energy, ψ = dmrg(H, ψ₀; nsweeps, maxdim, cutoff)

    updensity = fill(0.0, N)
    dndensity = fill(0.0, N)
    for i in 1:N
        orthogonalize!(ψ, i)
        ψᵢ = dag(prime(ψ[i], "Site"))
        updensity[i] = scalar(ψᵢ * op(sites, "Nup", i) * ψ[i])
        dndensity[i] = scalar(ψᵢ * op(sites, "Ndn", i) * ψ[i])
    end

    println("Up Density:")
    for ind in 1:N
        println("$ind $(updensity[ind])")
    end
    println()

    println("Dn Density")
    for ind in 1:N
        println("$ind $(dndensity[ind])")
    end
    println()

    println("Total Density")
    for ind in 1:N
        println("$ind $(updensity[ind] + dndensity[ind])")
    end
    println()

    println("\nGround State Energy = $energy")
end