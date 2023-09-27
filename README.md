# AD_TRACK
The transfer maps of Drift, Quad, Bend are in transfer_map_zygote.jl, corresponding to DriftPass, QuadLinearPass, BendLinearPass in AT.
The transfer maps are calculated using Zygote.Buffer that allows array mutation.

Use twissline to calculate optics of the lattice.
Use linepass1 to track multiple particles.
