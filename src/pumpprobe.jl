using DimensionalData
using Statistics

export delaytime_fs, position_mm_to_delaytime_fs
export tarpes_evolution

delaytime_fs(mirror_move_μm::Real) = 3.335640951981521 * mirror_move_μm
position_mm_to_delaytime_fs(position_mm::Real) = delaytime_fs(2 * position_mm * 1e3)

"""
    tarpes_evolution(A, delay_time=0.0, evolution_at=0.0; stack_dim=:delay, vertical_dim=:eV)
    tarpes_evolution(A, delay_index::Integer, evolution_at=0.0; ...)

Return a snapshot and the temporal evolution extracted from a 3D time-resolved ARPES dataset.

Arguments
- A::AbstractDimArray{T,3}: 3D tr-ARPES data with named dimensions (DimensionalData).
- delay_time::Real: Delay time (same units as the `stack_dim`) for the snapshot. Alternatively use delay_index.
- evolution_at::Real or Tuple{center, width}: If a scalar, select the nearest coordinate along the non-dispersion axis. If a (center, width) tuple, average over the window centered at `center` with the given `width`.
- stack_dim::Symbol: Name of the time/delay dimension (default :delay).
- vertical_dim::Symbol: Name of the vertical (energy) dimension (default :eV).
- full_temporal::Bool: If true, include the full time range for the temporal evolution; otherwise include only times <= `delay_time`.

Returns
- arpes_data_at_delay::AbstractDimArray{T,2}: 2D slice (non-dispersion × vertical) at the chosen delay.
- temporal_evolution_data::AbstractDimArray{T,2}: 2D array (time × vertical) with intensity evolution at the selected non-dispersion position(s).

Notes
- Selection uses DimensionalData indexing (Dim{...}(Near/Between)).
"""
function tarpes_evolution(
    A::AbstractDimArray{T,3} where {T},
    delay_time::Real = 0.0,
    evolution_at::Union{Real,Tuple{<:Real,<:Real}} = 0.0;
    stack_dim::Symbol = :delay,
    vertical_dim::Symbol = :eV,
    full_temporal::Bool = false,
)

    non_dispersion_axis = otherdims(A, (vertical_dim, stack_dim)) |> first |> name
    arpes_data_at_delay = A[Dim{stack_dim}(Near(delay_time))]

    temporal_evolution_data =
        _build_slice_data(A, non_dispersion_axis, evolution_at) |>
        x->permutedims(x, (stack_dim, vertical_dim)) |>
           x -> x[Dim{stack_dim}(Between(-Inf, full_temporal ? Inf : delay_time))]

    return arpes_data_at_delay, temporal_evolution_data
end

function tarpes_evolution(
    A::AbstractDimArray{T,3} where {T},
    delay_index::Integer,
    evolution_at::Union{Real,Tuple{<:Real,<:Real}} = 0.0;
    stack_dim::Symbol = :delay,
    vertical_dim::Symbol = :eV,
    full_temporal::Bool = false,
)
    delay_time = dims(A, stack_dim)[delay_index] |> float
    return tarpes_evolution(
        A,
        delay_time,
        evolution_at;
        stack_dim = stack_dim,
        vertical_dim = vertical_dim,
        full_temporal = full_temporal,
    )
end

"""
    _build_slice_data(A, non_dispersion_axis, evolution_at::Real)
    _build_slice_data(A, non_dispersion_axis, evolution_at::Tuple{center,width})

Internal helper. For a scalar `evolution_at`, select the nearest coordinate along
`non_dispersion_axis`. For a `(center, width)` tuple, average A over the window
[center - width/2, center + width/2] along `non_dispersion_axis` and drop that dimension.
"""
function _build_slice_data(
    A::AbstractDimArray{T,3} where {T},
    non_dispersion_axis::Symbol,
    evolution_at::Real,
)
    return A[Dim{non_dispersion_axis}(Near(evolution_at))]
end

function _build_slice_data(
    A::AbstractDimArray{T,3} where {T},
    non_dispersion_axis::Symbol,
    evolution_at::Tuple{<:Real,<:Real},
)
    evolution_left = evolution_at[1] - evolution_at[2]/2
    evolution_right = evolution_at[1] + evolution_at[2]/2

    sliced =
        A[Dim{non_dispersion_axis}(Between(evolution_left, evolution_right))] |>
        x ->
            mean(x, dims = non_dispersion_axis) |>
            x -> dropdims(x, dims = non_dispersion_axis)

    return sliced
end
