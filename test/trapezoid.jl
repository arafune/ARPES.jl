using Test
using ARPES
using ARPES: Z
using DimensionalData
using DimensionalData: At, dims, lookup, metadata, rebuild

if !@isdefined(test_ARPESData)
    include("fixture.jl")
end

const TEST_TRAPEZOID_CORNERS = (
    (eV = 6.0, phi = -5.0),
    (eV = 6.0, phi = 9.0),
    (eV = 8.0, phi = -3.0),
    (eV = 8.0, phi = 7.0),
)
const TEST_BASE_CORNERS = (-5.0, 9.0)

function trapezoid_test_slice()
    return test_ARPESData()[Z(At(1))]
end

function trapezoid_edges(E, corners)
    UL, UR, LL, LR = corners
    slope_left = (UL.phi - LL.phi) / (UL.eV - LL.eV)
    slope_right = (UR.phi - LR.phi) / (UR.eV - LR.eV)
    left = slope_left * (E - UL.eV) + UL.phi
    right = slope_right * (E - UR.eV) + UR.phi
    return left, right
end

function rectangle_to_trapezoid_phi(p, E, base_corners, corners)
    left, right = trapezoid_edges(E, corners)
    left_rect, right_rect = extrema(base_corners)
    c = (p - left_rect) / (right_rect - left_rect)
    return left + c * (right - left)
end

function trapezoid_to_rectangle_phi(p, E, corners, base_corners)
    left, right = trapezoid_edges(E, corners)
    left_rect, right_rect = extrema(base_corners)
    c = (p - left) / (right - left)
    return left_rect + c * (right_rect - left_rect)
end

function rectangle_reference_data(A)
    phi_vals = collect(lookup(A, :phi))
    eV_vals = collect(lookup(A, :eV))
    return [phi + 10eV for phi in phi_vals, eV in eV_vals]
end

function trapezoid_reference_data(A, corners, base_corners)
    phi_vals = collect(lookup(A, :phi))
    eV_vals = collect(lookup(A, :eV))
    out = Matrix{Float64}(undef, length(phi_vals), length(eV_vals))

    for (j, E) in pairs(eV_vals), (i, p) in pairs(phi_vals)
        left, right = trapezoid_edges(E, corners)
        if left <= p <= right
            rect_phi = trapezoid_to_rectangle_phi(p, E, corners, base_corners)
            out[i, j] = rect_phi + 10E
        else
            out[i, j] = NaN
        end
    end

    return out
end

function trapezoid_full_reference_data(A, corners, base_corners)
    phi_vals = collect(lookup(A, :phi))
    eV_vals = collect(lookup(A, :eV))
    return [
        trapezoid_to_rectangle_phi(phi, eV, corners, base_corners) + 10eV for
        phi in phi_vals, eV in eV_vals
    ]
end

function test_matrix_equal_or_nan(actual, expected; atol = 1e-8)
    @test size(actual) == size(expected)
    for idx in eachindex(actual, expected)
        if isnan(expected[idx])
            @test isnan(actual[idx])
        else
            @test isapprox(actual[idx], expected[idx]; atol = atol)
        end
    end
end

@testset "_phi_to_phi maps between rectangle and trapezoid" begin
    rect_phi = [-5.0, 2.0, 9.0]
    Ek = [6.0, 7.0, 8.0]
    trap_phi = [
        rectangle_to_trapezoid_phi(rect_phi[i], Ek[i], TEST_BASE_CORNERS, TEST_TRAPEZOID_CORNERS) for
        i in eachindex(rect_phi)
    ]

    @test ARPES._phi_to_phi(rect_phi, Ek, TEST_BASE_CORNERS, TEST_TRAPEZOID_CORNERS) ≈ trap_phi
    @test ARPES._phi_to_phi(
        trap_phi,
        Ek,
        TEST_TRAPEZOID_CORNERS,
        TEST_BASE_CORNERS,
    ) ≈ rect_phi
end

@testset "trapezoid converts rectangle data to trapezoid coordinates" begin
    A = trapezoid_test_slice()
    rect_data = rectangle_reference_data(A)
    A_rect = rebuild(A; data = rect_data)

    converted = trapezoid(A_rect, TEST_BASE_CORNERS, TEST_TRAPEZOID_CORNERS)
    expected = trapezoid_reference_data(A, TEST_TRAPEZOID_CORNERS, TEST_BASE_CORNERS)

    @test collect(lookup(converted, :phi)) == collect(lookup(A_rect, :phi))
    @test collect(lookup(converted, :eV)) == collect(lookup(A_rect, :eV))
    @test metadata(converted) == metadata(A_rect)
    @test metadata(dims(converted, :phi)) == metadata(dims(A_rect, :phi))
    @test metadata(dims(converted, :eV)) == metadata(dims(A_rect, :eV))
    test_matrix_equal_or_nan(parent(converted), expected)
end

@testset "trapezoid converts trapezoid data back to rectangle coordinates" begin
    A = trapezoid_test_slice()
    input_phi_vals = collect(lookup(A, :phi))
    eV_vals = collect(lookup(A, :eV))
    step_phi = abs(input_phi_vals[2] - input_phi_vals[1])
    trap_data = trapezoid_full_reference_data(A, TEST_TRAPEZOID_CORNERS, TEST_BASE_CORNERS)
    A_trapezoid = rebuild(A; data = trap_data)

    converted = trapezoid(A_trapezoid, TEST_TRAPEZOID_CORNERS, TEST_BASE_CORNERS)
    expected_phi = collect(
        range(
            start = minimum(
                trapezoid_to_rectangle_phi(
                    minimum(input_phi_vals),
                    E,
                    TEST_TRAPEZOID_CORNERS,
                    TEST_BASE_CORNERS,
                ) for E in eV_vals
            ),
            step = step_phi,
            stop = maximum(
                trapezoid_to_rectangle_phi(
                    maximum(input_phi_vals),
                    E,
                    TEST_TRAPEZOID_CORNERS,
                    TEST_BASE_CORNERS,
                ) for E in eV_vals
            ),
        ),
    )
    expected = Matrix{Float64}(undef, length(expected_phi), length(eV_vals))
    min_input_phi = minimum(input_phi_vals)
    max_input_phi = maximum(input_phi_vals)
    for (j, E) in pairs(eV_vals), (i, phi) in pairs(expected_phi)
        trap_phi = rectangle_to_trapezoid_phi(phi, E, TEST_BASE_CORNERS, TEST_TRAPEZOID_CORNERS)
        if min_input_phi <= trap_phi <= max_input_phi
            expected[i, j] = phi + 10E
        else
            expected[i, j] = NaN
        end
    end

    @test collect(lookup(converted, :phi)) ≈ expected_phi
    @test collect(lookup(converted, :eV)) == eV_vals
    @test metadata(converted) == metadata(A_trapezoid)
    test_matrix_equal_or_nan(parent(converted), expected)
end

@testset "trapezoid validates input geometry" begin
    A = trapezoid_test_slice()
    rect_data = rectangle_reference_data(A)
    A_rect = rebuild(A; data = rect_data)

    phi_vals = collect(lookup(A_rect, :phi))
    irregular_phi = copy(phi_vals)
    irregular_phi[2] += 0.1
    irregular_dim = rebuild(dims(A_rect, :phi); val = irregular_phi)
    irregular_A = rebuild(A_rect; dims = (irregular_dim, dims(A_rect, :eV)))

    bad_corners = (
        (eV = 6.0, phi = -5.0),
        (eV = 6.5, phi = 9.0),
        (eV = 8.0, phi = -3.0),
        (eV = 8.0, phi = 7.0),
    )

    @test_throws AssertionError trapezoid(irregular_A, TEST_BASE_CORNERS, TEST_TRAPEZOID_CORNERS)
    @test_throws AssertionError trapezoid(A_rect, (9.0, -5.0), TEST_TRAPEZOID_CORNERS)
    @test_throws AssertionError trapezoid(A_rect, TEST_BASE_CORNERS, bad_corners)
end
