using DimensionalData.Dimensions: @dim

#--- Dimension names used for ARPES ---
@dim kx "kx  ( Å⁻¹ )"
@dim ky "ky  ( Å⁻¹ )"
@dim kz "kz  ( Å⁻¹ )"
@dim psi "ψ  ( degrees )"
@dim phi "Φ  ( degrees )"
@dim eV "eV"
@dim detector_ch
@dim ch2
@dim delay
@dim spin

"""
type representing the unit of intensity in ARPES data.
"""
abstract type IntensityUnit end
struct Counts <: IntensityUnit end
struct CPS <: IntensityUnit end

"""
type representing the analyzer configuration defined in Rev. Sci. Instrum. **89**, 043903 (2018).
"""
abstract type AnalyzerConfiguration end
struct TypeI <: AnalyzerConfiguration end
struct TypeII <: AnalyzerConfiguration end
struct TypeIp <: AnalyzerConfiguration end
struct TypeIIp <: AnalyzerConfiguration end

"""
    @enum EnergyDefinition

Enumeration of possible energy definitions used in ARPES data analysis.

- `BindingEnergy`: Electron binding energy relative to the Fermi level.
- `FinalStateEnergy`: Energy of the electron in the final state after photoemission (Referred to the Fermi level).
- `KineticEnergy`: Kinetic energy of the emitted electron (Referred to the vacuum level of the sample).
- `IntermediateEnergy`: Energy in an intermediate state (e.g., in pump-probe experiments).
"""
@enum EnergyDefinition begin
    BindingEnergy
    FinalStateEnergy
    KineticEnergy
    IntermediateEnergy
end
