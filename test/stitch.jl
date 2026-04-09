using Test
using DimensionalData
using ARPES


function build_arrays()
    x = range(1.0, stop = 8.0, length = 20)
    y = collect(range(5, stop = 9, length = 10))
    data = collect(1.0:200) |> d->reshape(d, 20, 10)
    A = DimArray(data, (X = x, Y = y))

    x = range(3.0, stop = 15.0, length = 20)
    y = collect(range(5, stop = 9, length = 10))
    data = collect(1.0:200) |> d->reshape(d, 10, 20)
    B = DimArray(data', (X = x, Y = y))

    x = range(13.0, stop = 15.0, length = 10)
    y = collect(range(5, stop = 9, length = 10))
    data = collect(1.0:100) |> d->reshape(d, 10, 10)
    C = DimArray(data', (X = x, Y = y))

    return A, B, C

end

function build_arrays3D()
    x = range(1.0, stop = 8.0, length = 20)
    y = collect(range(5, stop = 9, length = 10))
    z = collect(range(1, stop = 5, length = 5))
    data = collect(1.0:1000) |> d->reshape(d, 20, 10, 5)
    A = DimArray(data, (X = x, Y = y, Z = z))

    x = range(7.0, stop = 15.0, length = 20)
    y = collect(range(5, stop = 9, length = 10))
    z = collect(range(1, stop = 5, length = 5))
    data = collect(1.0:1000) * 2.1 |> d->reshape(d, 20, 10, 5)
    B = DimArray(data, (X = x, Y = y, Z = z))
    return A, B
end

