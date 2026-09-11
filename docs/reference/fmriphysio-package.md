# fmriphysio: Fast Design-Free Physiological Denoising for Functional MRI

Removes physiological noise from BOLD fMRI time series without a task
design. A fast temporal-difference energy map
([`compute_wNN_diff_energy()`](compute_wNN_diff_energy.md)) separates
likely neuronal voxels from non-neuronal voxels (vasculature,
cerebrospinal fluid). Temporally predictable components are extracted
from a low-rank subspace, retained when they load more strongly on
non-neuronal than on neuronal voxels, and projected out of the data
([`fast_phy_denoise()`](fast_phy_denoise.md)). A CompCor baseline
([`compcor_denoise()`](compcor_denoise.md)) and design-free
quality-control metrics
([`compute_design_free_qc()`](compute_design_free_qc.md)) are provided
for comparison.

## References

Churchill, N. W., & Strother, S. C. (2013). PHYCAA+: An optimized,
adaptive procedure for measuring and controlling physiological noise in
BOLD fMRI. *NeuroImage*, 82, 306–325.
[doi:10.1016/j.neuroimage.2013.05.102](https://doi.org/10.1016/j.neuroimage.2013.05.102)

Churchill, N. W., Yourganov, G., Spring, R., Rasmussen, P. M., Lee, W.,
Ween, J. E., & Strother, S. C. (2012). PHYCAA: Data-driven measurement
and removal of physiological noise in BOLD fMRI. *NeuroImage*, 59(2),
1299–1314.
[doi:10.1016/j.neuroimage.2011.08.021](https://doi.org/10.1016/j.neuroimage.2011.08.021)

Behzadi, Y., Restom, K., Liau, J., & Liu, T. T. (2007). A component
based noise correction method (CompCor) for BOLD and perfusion based
fMRI. *NeuroImage*, 37(1), 90–101.
[doi:10.1016/j.neuroimage.2007.04.042](https://doi.org/10.1016/j.neuroimage.2007.04.042)

Dong, Y., & Qin, S. J. (2018). A novel dynamic PCA algorithm for dynamic
data modeling and process monitoring. *Journal of Process Control*, 67,
1–11.
[doi:10.1016/j.jprocont.2017.05.002](https://doi.org/10.1016/j.jprocont.2017.05.002)

## See also

Useful links:

- <https://github.com/bbuchsbaum/fmriphysio>

- Report bugs at <https://github.com/bbuchsbaum/fmriphysio/issues>

## Author

**Maintainer**: Bradley Buchsbaum <brad.buchsbaum@gmail.com>
