"""
    Format

A module for handling different file formats.

This module provides functionality for converting measurement data into a DimArray.
It includes a structure for organizing code related to the various file formats that
may be used in the application.

Each format has its own implementation and functions for reading data.
Data is read in its original form as much as possible.
Beamline-specific rules should be described in Location modules under the `location` directory.

It includes the `itx.jl` file, which contains the implementation for the ITX format.
This module can be extended to include additional formats as needed.mo

In the `Format` module, functions for reading data from specific formats should be defined,
and any necessary conversions to the standard `DimArray` format should be implemented.
The module serves as a central place for handling different file formats and their associated data processing logic.
"""
module Format
include("itx.jl")
end
