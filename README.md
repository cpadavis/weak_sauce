# Weak Sensor Anomalies Under Careful Examination (weak_sauce)

## Summary

This package, described in detail in [Baumer, et al. 2017](), provides methods to fit models of distorted pixel grids to flat-field data, and compute their corresponding systematic impacts. For examples of how the code can be used, see the included [demo notebook]().

## Installation

After cloning, there are two submodules to compile before the code is ready to use.

1. Within `code/weak_sauce/r2d/`, compile the pixel flux integration module with:
> $> make

2. Within `code/weak_sauce/adaptive_moments/`, compile the adaptive moments module with:
> $> python adaptive_moments_setup.py build_ext --inplace
