The project is about studying an incompressible fluid flow near a porous material.

The idea:
- Divide the domain in 3 regions:
    1. Upper of undisturbed fluid flow (DNS)
    1. Middle with porous material solving the inner fluid flow (DNS)
    1. Lower with porous material to impose fisical BCs, no need to  study the flow in detail (Volume Avaraged N-S eqs aka VANS)

- Implement a solver for each region:
    1. Fourier-Fourier-Chebichev semi-spectral // pseudo-spectral method here, some tests for efficency have to be made.  
    1. Finite differences for the momentum eq and FFT approach for the laplace eq for the pressure. (Here an Immersed Buondary Method needs to be implemented or at least is the most "easy" way I saw people treat this problem in the literature, I have made some  tests alreay on simple 1D case so the key ideas should be there (hopefully))
    1. Needs to be studied better the implementation here.

- Go parallel (and very fast possibly)
     1. This first domain shoud be a piece of cake to make parallel (GPU acceleration should  be great)
     1. This second layer has to be treated with some care due to the IBM approach, but ideally should not be impossible at least
     1. Again needs to be understood the approach here, but ideally in the worst case scenario is a linear  system to solve so should not be that bad.

- What we have: 
    1. An already functioning DNS code that I contributed to produce in the context of the HPSC for Aerospace course (help by profs Auteri and Piscaglia for the HPC engineering master degree), should be a good startig point at least: 
        - the tests I made on my machine suggest speed in the order of 10^-8 seconds/(N_nodes * N_processors * Timestep)
        - the code uses already FD and FFTs as suggested above, is fully explicit and uses RK3 for time integration 
        - the code is implemented for parallel run on CPU with MPI, but I'd like to rewrite it and to make it GPU-friendly
    1. Lot of passion and will to make mistakes and do a first code of this type to make of public use as the study of porous material is becoming more and more relevant in Aerospace industry due to theyr ability of enabling turbulence transition at lower Reynolds numbers.
    1. Prof Cortelezzi has a very good experience on porous materials so I'm confident that his guide will be helpful in the validation of the code.
