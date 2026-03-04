"""Bilinear interpolation."""
using Base.Threads

"""
    fast_bilinear_interpolate(
        desired_pos_dim0::AbstractArray,
        desired_pos_dim1::AbstractArray,
        orig_coords_dim0::AbstractVector,
        orig_coords_dim1::AbstractVector,
        orig_values::AbstractMatrix,
    ) -> AbstractArray

Performs bilinear interpolation on a 2D grid for arbitrary (not necessarily rectilinear) axes.

# Arguments
- `desired_pos_dim0::AbstractArray`: Desired positions along the first dimension (x-axis) to interpolate.
- `desired_pos_dim1::AbstractArray`: Desired positions along the second dimension (y-axis) to interpolate.
- `orig_coords_dim0::AbstractVector`: Original grid coordinates along the first dimension (x-axis), must be sorted.
- `orig_coords_dim1::AbstractVector`: Original grid coordinates along the second dimension (y-axis), must be sorted.
- `orig_values::AbstractMatrix`: Matrix of values defined on the original grid.

# Returns
- `AbstractArray`: Interpolated values at the desired positions, with the same shape as the input position arrays.

Returns `NaN` for points outside the original grid.

Threaded for performance.
"""
using Base.Threads

function fast_bilinear_interpolate(
    desired_pos_dim0::AbstractArray,
    desired_pos_dim1::AbstractArray,
    orig_coords_dim0::AbstractVector,
    orig_coords_dim1::AbstractVector,
    orig_values::AbstractMatrix,
)
    # Check axis order and set flags
    rev_dim0 = orig_coords_dim0[1] > orig_coords_dim0[end]
    rev_dim1 = orig_coords_dim1[1] > orig_coords_dim1[end]

    coords_dim0 = rev_dim0 ? reverse(orig_coords_dim0) : orig_coords_dim0
    coords_dim1 = rev_dim1 ? reverse(orig_coords_dim1) : orig_coords_dim1
    values = orig_values
    if rev_dim0
        values = reverse(values, dims = 1)
    end
    if rev_dim1
        values = reverse(values, dims = 2)
    end

    len_dim0 = length(coords_dim0)
    len_dim1 = length(coords_dim1)
    desired_shape = size(desired_pos_dim0)
    xs = vec(desired_pos_dim0)
    ys = vec(desired_pos_dim1)
    n_points = length(xs)
    result = Vector{eltype(orig_values)}(undef, n_points)

    @threads for idx = 1:n_points
        x = xs[idx]
        y = ys[idx]

        x1_idx = searchsortedfirst(coords_dim0, x) - 1
        x2_idx = x1_idx + 1
        y1_idx = searchsortedfirst(coords_dim1, y) - 1
        y2_idx = y1_idx + 1

        if x1_idx < 1 || x2_idx > len_dim0 || y1_idx < 1 || y2_idx > len_dim1
            result[idx] = NaN
            continue
        end

        @inbounds begin
            x1 = coords_dim0[x1_idx]
            x2 = coords_dim0[x2_idx]
            y1 = coords_dim1[y1_idx]
            y2 = coords_dim1[y2_idx]

            Q11 = values[x1_idx, y1_idx]
            Q12 = values[x1_idx, y2_idx]
            Q21 = values[x2_idx, y1_idx]
            Q22 = values[x2_idx, y2_idx]

            denom = (x2 - x1) * (y2 - y1)
            result[idx] =
                (
                    Q11 * (x2 - x) * (y2 - y) +
                    Q21 * (x - x1) * (y2 - y) +
                    Q12 * (x2 - x) * (y - y1) +
                    Q22 * (x - x1) * (y - y1)
                ) / denom
        end
    end

    return reshape(result, desired_shape)
end


"""
    fast_bilinear_interpolate_rectilinear(
        desired_pos_dim0::AbstractArray,
        desired_pos_dim1::AbstractArray,
        orig_coords_dim0::AbstractVector,
        orig_coords_dim1::AbstractVector,
        orig_values::AbstractMatrix,
    ) -> AbstractArray

Performs fast bilinear interpolation on a 2D rectilinear grid (uniformly spaced axes).

# Arguments
- `desired_pos_dim0::AbstractArray`: Desired positions along the first dimension (x-axis) to interpolate.
- `desired_pos_dim1::AbstractArray`: Desired positions along the second dimension (y-axis) to interpolate.
- `orig_coords_dim0::AbstractVector`: Original grid coordinates along the first dimension (x-axis), must be sorted and uniformly spaced.
- `orig_coords_dim1::AbstractVector`: Original grid coordinates along the second dimension (y-axis), must be sorted and uniformly spaced.
- `orig_values::AbstractMatrix`: Matrix of values defined on the original grid.

# Returns
- `AbstractArray`: Interpolated values at the desired positions, with the same shape as the input position arrays.

Returns `NaN` for points outside the original grid.

Threaded for performance.
"""
function fast_bilinear_interpolate_rectilinear(
    desired_pos_dim0::AbstractArray,
    desired_pos_dim1::AbstractArray,
    orig_coords_dim0::AbstractVector,
    orig_coords_dim1::AbstractVector,
    orig_values::AbstractMatrix,
)
    rev_dim0 = orig_coords_dim0[1] > orig_coords_dim0[end]
    rev_dim1 = orig_coords_dim1[1] > orig_coords_dim1[end]

    coords_dim0 = rev_dim0 ? reverse(orig_coords_dim0) : orig_coords_dim0
    coords_dim1 = rev_dim1 ? reverse(orig_coords_dim1) : orig_coords_dim1
    values = orig_values
    if rev_dim0
        values = reverse(values, dims = 1)
    end
    if rev_dim1
        values = reverse(values, dims = 2)
    end

    len_dim0 = length(coords_dim0)
    len_dim1 = length(coords_dim1)
    desired_shape = size(desired_pos_dim0)
    xs = vec(desired_pos_dim0)
    ys = vec(desired_pos_dim1)
    n_points = length(xs)
    result = Vector{eltype(orig_values)}(undef, n_points)

    step_dim0 = (coords_dim0[end] - coords_dim0[1]) / (len_dim0 - 1)
    step_dim1 = (coords_dim1[end] - coords_dim1[1]) / (len_dim1 - 1)
    x0 = coords_dim0[1]
    y0 = coords_dim1[1]

    @threads for idx = 1:n_points
        x = xs[idx]
        y = ys[idx]

        x1_idx = floor(Int, (x - x0) / step_dim0) + 1
        x2_idx = x1_idx + 1
        y1_idx = floor(Int, (y - y0) / step_dim1) + 1
        y2_idx = y1_idx + 1

        if x1_idx < 1 || x2_idx > len_dim0 || y1_idx < 1 || y2_idx > len_dim1
            result[idx] = NaN
            continue
        end

        @inbounds begin
            x1 = coords_dim0[x1_idx]
            x2 = coords_dim0[x2_idx]
            y1 = coords_dim1[y1_idx]
            y2 = coords_dim1[y2_idx]

            Q11 = values[x1_idx, y1_idx]
            Q12 = values[x1_idx, y2_idx]
            Q21 = values[x2_idx, y1_idx]
            Q22 = values[x2_idx, y2_idx]

            denom = (x2 - x1) * (y2 - y1)
            result[idx] =
                (
                    Q11 * (x2 - x) * (y2 - y) +
                    Q21 * (x - x1) * (y2 - y) +
                    Q12 * (x2 - x) * (y - y1) +
                    Q22 * (x - x1) * (y - y1)
                ) / denom
        end
    end

    return reshape(result, desired_shape)
end

function fast_trilinear_interpolate(
    desired_pos_dim0::AbstractArray,
    desired_pos_dim1::AbstractArray,
    desired_pos_dim2::AbstractArray,
    orig_coords_dim0::AbstractVector,
    orig_coords_dim1::AbstractVector,
    orig_coords_dim2::AbstractVector,
    orig_values::AbstractArray{T,3};
    progress_proxy = nothing,
) where {T}
    # 軸の向きを確認し、必要に応じて反転（コピーを避けるため後でViewを使用検討）
    rev0 = orig_coords_dim0[1] > orig_coords_dim0[end]
    rev1 = orig_coords_dim1[1] > orig_coords_dim1[end]
    rev2 = orig_coords_dim2[1] > orig_coords_dim2[end]

    coords0 = rev0 ? reverse(orig_coords_dim0) : orig_coords_dim0
    coords1 = rev1 ? reverse(orig_coords_dim1) : orig_coords_dim1
    coords2 = rev2 ? reverse(orig_coords_dim2) : orig_coords_dim2

    # 値の配列を座標系に合わせて反転
    values = orig_values
    if rev0
        ;
        values = reverse(values, dims = 1);
    end
    if rev1
        ;
        values = reverse(values, dims = 2);
    end
    if rev2
        ;
        values = reverse(values, dims = 3);
    end

    len0, len1, len2 = length(coords0), length(coords1), length(coords2)
    desired_shape = size(desired_pos_dim0)

    # 計算用にフラット化
    xs = vec(desired_pos_dim0)
    ys = vec(desired_pos_dim1)
    zs = vec(desired_pos_dim2)

    n_points = length(xs)
    result = Vector{T}(undef, n_points)

    @threads for idx = 1:n_points
        # 進捗更新（必要な場合）
        if !isnothing(progress_proxy) && (idx % 100 == 0 || idx == n_points)
            # 注: JuliaのThreads環境下でのProgress更新は
            # ライブラリ側のスレッドセーフ性に依存します
        end

        x, y, z = xs[idx], ys[idx], zs[idx]

        # 1ベースインデックスへの変換を考慮した探索
        i1 = searchsortedfirst(coords0, x) - 1
        i2 = i1 + 1
        j1 = searchsortedfirst(coords1, y) - 1
        j2 = j1 + 1
        k1 = searchsortedfirst(coords2, z) - 1
        k2 = k1 + 1

        # 境界チェック
        if i1 < 1 || i2 > len0 || j1 < 1 || j2 > len1 || k1 < 1 || k2 > len2
            result[idx] = NaN
            continue
        end

        @inbounds begin
            x1, x2 = coords0[i1], coords0[i2]
            y1, y2 = coords1[j1], coords1[j2]
            z1, z2 = coords2[k1], coords2[k2]

            # 8つの隣接点の値を抽出
            q111 = values[i1, j1, k1]
            q112 = values[i1, j1, k2]
            q121 = values[i1, j2, k1]
            q122 = values[i1, j2, k2]
            q211 = values[i2, j1, k1]
            q212 = values[i2, j1, k2]
            q221 = values[i2, j2, k1]
            q222 = values[i2, j2, k2]

            # 三線形補間
            # 分母の計算
            denom = (x2 - x1) * (y2 - y1) * (z2 - z1)

            val =
                (
                    q111 * (x2 - x) * (y2 - y) * (z2 - z) +
                    q211 * (x - x1) * (y2 - y) * (z2 - z) +
                    q121 * (x2 - x) * (y - y1) * (z2 - z) +
                    q221 * (x - x1) * (y - y1) * (z2 - z) +
                    q112 * (x2 - x) * (y2 - y) * (z - z1) +
                    q212 * (x - x1) * (y2 - y) * (z - z1) +
                    q122 * (x2 - x) * (y - y1) * (z - z1) +
                    q222 * (x - x1) * (y - y1) * (z - z1)
                ) / denom

            result[idx] = val
        end
    end

    return reshape(result, desired_shape)
end

function fast_trilinear_interpolate_rectilinear(
    desired_pos_dim0::AbstractArray,
    desired_pos_dim1::AbstractArray,
    desired_pos_dim2::AbstractArray,
    orig_coords_dim0::AbstractVector,
    orig_coords_dim1::AbstractVector,
    orig_coords_dim2::AbstractVector,
    orig_values::AbstractArray{T,3};
    progress_proxy = nothing,
) where {T}
    # 軸の反転処理
    rev0 = orig_coords_dim0[1] > orig_coords_dim0[end]
    rev1 = orig_coords_dim1[1] > orig_coords_dim1[end]
    rev2 = orig_coords_dim2[1] > orig_coords_dim2[end]

    coords0 = rev0 ? reverse(orig_coords_dim0) : orig_coords_dim0
    coords1 = rev1 ? reverse(orig_coords_dim1) : orig_coords_dim1
    coords2 = rev2 ? reverse(orig_coords_dim2) : orig_coords_dim2

    values = orig_values
    if rev0
        ;
        values = reverse(values, dims = 1);
    end
    if rev1
        ;
        values = reverse(values, dims = 2);
    end
    if rev2
        ;
        values = reverse(values, dims = 3);
    end

    len0, len1, len2 = length(coords0), length(coords1), length(coords2)
    desired_shape = size(desired_pos_dim0)

    xs = vec(desired_pos_dim0)
    ys = vec(desired_pos_dim1)
    zs = vec(desired_pos_dim2)
    n_points = length(xs)
    result = Vector{T}(undef, n_points)

    # ステップサイズの計算
    step0 = (coords0[end] - coords0[1]) / (len0 - 1)
    step1 = (coords1[end] - coords1[1]) / (len1 - 1)
    step2 = (coords2[end] - coords2[1]) / (len2 - 1)

    @threads for idx = 1:n_points
        x, y, z = xs[idx], ys[idx], zs[idx]

        # インデックスの計算 (Pythonの int((x - start) / step) に相当)
        # Juliaでは1ベースなので、最後に +1 する
        i1 = floor(Int, (x - coords0[1]) / step0) + 1
        i2 = i1 + 1
        j1 = floor(Int, (y - coords1[1]) / step1) + 1
        j2 = j1 + 1
        k1 = floor(Int, (z - coords2[1]) / step2) + 1
        k2 = k1 + 1

        # 境界チェック
        if i1 < 1 || i2 > len0 || j1 < 1 || j2 > len1 || k1 < 1 || k2 > len2
            result[idx] = NaN
            continue
        end

        @inbounds begin
            x1, x2 = coords0[i1], coords0[i2]
            y1, y2 = coords1[j1], coords1[j2]
            z1, z2 = coords2[k1], coords2[k2]

            q111 = values[i1, j1, k1]
            q112 = values[i1, j1, k2]
            q121 = values[i1, j2, k1]
            q122 = values[i1, j2, k2]
            q211 = values[i2, j1, k1]
            q212 = values[i2, j1, k2]
            q221 = values[i2, j2, k1]
            q222 = values[i2, j2, k2]

            denom = (x2 - x1) * (y2 - y1) * (z2 - z1)

            val =
                (
                    q111 * (x2 - x) * (y2 - y) * (z2 - z) +
                    q211 * (x - x1) * (y2 - y) * (z2 - z) +
                    q121 * (x2 - x) * (y - y1) * (z2 - z) +
                    q221 * (x - x1) * (y - y1) * (z2 - z) +
                    q112 * (x2 - x) * (y2 - y) * (z - z1) +
                    q212 * (x - x1) * (y2 - y) * (z - z1) +
                    q122 * (x2 - x) * (y - y1) * (z - z1) +
                    q222 * (x - x1) * (y - y1) * (z - z1)
                ) / denom

            result[idx] = val
        end
    end

    return reshape(result, desired_shape)
end
