### Example code for the Hubbard model on a tri-angular lattice
using ITensors
let
    Ny = 6
    Nx = 12
    N = Nx * Ny
    
    sites = siteinds("S=1/2", N; conserve_qns = true)

    # Obtain an array of LatticeBond structs
    lattice = triangular_lattice(Nx, Ny, yperiodic = false)

    
    
    # Define the Heisenberg spin Hamiltonian on this lattice
    # Implement the Hubbard model [To Do].
    ampo = OpSum()
    for bond in lattice
        ampo .+= 0.5, "S+", bond.s1, "S-", bond.s2
        ampo .+= 0.5, "S-", bond.s1, "S+", bond.s2
        ampo .+=      "Sz", bond.s1, "Sz", bond.s2
    end
    H = MPO(ampo, sites)

    state = [isodd(n) ? "Up" : "Dn" for n = 1:N]
    ψ₀ = randomMPS(sites, state, 20)

    sweeps = Sweeps(10)
    maxdim!(sweeps, 20, 70, 100, 100, 200, 400, 800)
    cutoff!(sweeps, 1E-8)
    @show sweeps

    energy, ψ = dmrg(H, ψ₀, sweeps)
    return
end