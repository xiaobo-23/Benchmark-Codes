### Example code for the Hubbard model on a tri-angular lattice
using ITensors
let
    Ny = 3
    Nx = 12
    N = Nx * Ny
    t = 1.0                 # The hopping amplitude in the Hubbard model
    U = 8.0                 # The on-site Coulomb interaction in the Hubbard model
    sites = siteinds("Electron", N; conserve_qns = true)

    # Obtain an array of LatticeBond structs
    lattice = triangular_lattice(Nx, Ny, yperiodic = true)

    
    # Define the Hubbard Hamiltonian on a triangular lattice
    ampo = OpSum()
    for bond in lattice
        ampo += -t, "Cdagup", bond.s1, "Cup", bond.s2
        ampo += -t, "Cdagup", bond.s2, "Cup", bond.s1
        ampo += -t, "Cdagdn", bond.s1, "Cdn", bond.s2
        ampo += -t, "Cdagdn", bond.s2, "Cdn", bond.s1
    end

    # On-site interaction
    for n in 1:N
        ampo += U, "Nupdn", N
    end

    H = MPO(ampo, sites)

    state = [isodd(n) ? "Up" : "Dn" for n = 1:N]
    ψ₀ = randomMPS(sites, state, 20)

    sweeps = Sweeps(20)
    maxdim!(sweeps, 20, 70, 100, 100, 200, 400, 800, 1000, 1500)
    cutoff!(sweeps, 1E-8)
    @show sweeps

    energy, ψ = dmrg(H, ψ₀, sweeps)
    @show t, U
    @show flux(ψ)
    @show maxlinkdim(ψ)
    @show energy
    return
end