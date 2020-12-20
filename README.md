# Weak Sensor Anomalies Under Careful Examination (weak_sauce)

## Summary

This package, described in detail in [Baumer, et al. 2017](https://arxiv.org/abs/1706.07400), provides methods to fit models of distorted pixel grids to flat-field data, and compute their corresponding systematic impacts. For examples of how the code can be used, see the included [demo notebook](https://github.com/cpadavis/weak_sauce/blob/master/notebooks/demo.ipynb).

## Installation

After cloning, there are two submodules to compile before the code is ready to use.

1. Within `code/weak_sauce/r3d/`, compile the pixel flux integration module with:
> $> make

2. Within `code/weak_sauce/adaptive_moments/`, compile the adaptive moments module with:
> $> python adaptive_moments_setup.py build_ext --inplace

3. Make sure you have [https://github.com/GalSim-developers/GalSim](galsim) installed (easiest done with `conda`).
> $> conda install galsim

4. You'll also need to install `sklearn, scipy, numpy, matplotlib`.

5. Make sure that the path to `code` is on your `$PYTHONPATH`, or alternatively in your code add the following lines:

> import sys
> sys.path.append("/path/to/weak_sauce/code")
