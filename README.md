# CanopyOptics.jl
<p align="center">
  
  <a href="https://RemoteSensingTools.github.io/CanopyOptics.jl/dev/">
    <img src="https://img.shields.io/badge/docs-latest-blue.svg"
         alt="Docs">
  </a>
  <a href="https://github.com/RemoteSensingTools/CanopyOptics.jl/blob/master/LICENSE">
    <img src="https://img.shields.io/github/license/RemoteSensingTools/CanopyOptics.jl"
         alt="License">
  </a>
  <a href="https://github.com/RemoteSensingTools/vSmartMOM.jl/commits/master">
    <img src="https://img.shields.io/github/commit-activity/y/RemoteSensingTools/CanopyOptics.jl"
         alt="Github Commit Frequency">
  </a>
</p>

Tools for the computation of canopy optical parameters, see docs for details (works in concert with vSmartMOM.jl)

Convention note: canopy Z matrices use `Z[i_out, j_in]` and normalize
BiLambertian leaves by `G(μ_in) * (R + T)`, matching the vSmartMOM
phase-matrix convention. Use `compute_Z_matrices(model, μ, LD, 0:m_max)` to
build a full Fourier stack.
