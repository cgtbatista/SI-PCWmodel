# Unraveling the Mechanical Behavior of Softwood Secondary Cell Walls through Atomistic Simulations

These files and scripts belong to the Supporting Information of the article: X-TITLE-X. AUTHORS Lucas N. Trentin, Amadeus C. S. de Alcântara, Carlos G T Batista, Munir S. Skaf ACS xxx https://doi.org/xxx.



### How to use the distance.jl script?

You may use or adapt the julia script to get the distance data of each one of the hemicellulose chains based on the cellulose macrofibril. On Julia, you should run something like:

```julia
includet("./distance.jl")

pdb, dcd = "/home/user/Documents/phd/sandbox/pcw.pdb", "/home/user/Documents/phd/sandbox/traj.dcd"
mindist(pdb, dcd, "XY12"; fstep=1, dist_cutoff=40.0)
```

The **pdb** and **dcd** are the PDB and the trajectory files. In this case, the function will walk over 1 frame of the trajectory and compute the minimum distances of XY12 xylan based on 40.0 Å cutoff. The XY12 xylan is the segname of the xylan. Since this process can take some time for a big system and trajectory, I suggest put all the mention code on a script (*e.g.* d_min.jl) and run Julia on multi-threads like `julia -t 8 d_min.jl`.