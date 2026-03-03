using Interpolations
using Base.Threads


function fast_bilinear_interpolate(
    desired_pos_dim0::AbstractArray,
    desired_pos_dim1::AbstractArray,
    orig_coords_dim0::AbstractVector,
    orig_coords_dim1::AbstractVector,
    orig_values::AbstractMatrix,
)
    # 1. 補間オブジェクト（Itp）の作成
    # Gridded(Linear()) は不等間隔なグリッドにも対応する線形補間を指定
    nodes = (orig_coords_dim0, orig_coords_dim1)
    itp = interpolate(nodes, orig_values, Gridded(Linear()))

    # 2. 範囲外を NaN に設定
    sitp = extrapolate(itp, NaN)

    # 3. 補間実行（ブロードキャストによる並列化または明示的なループ）
    # reshapeを自動で扱うため、入力と同じ形状で返せる
    return sitp.(desired_pos_dim0, desired_pos_dim1)
end


function fast_bilinear_interpolate_rectilinear(
    desired_pos_dim0,
    desired_pos_dim1,
    orig_coords_dim0,
    orig_coords_dim1,
    orig_values,
)
    # 1. 補間オブジェクトの作成
    # 等間隔（Regular Grid）の場合は BSpline(Linear()) が最適
    itp = interpolate(orig_values, BSpline(Linear()))

    # 2. 補間オブジェクトに実際の物理座標（軸）を割り当てる
    # これにより、任意の (x, y) でのアクセスが可能になる
    sitp = scale(itp, orig_coords_dim0, orig_coords_dim1)

    # 3. 範囲外の挙動を NaN に設定
    eitp = extrapolate(sitp, NaN)

    # 4. ブロードキャスト（.演算子）で計算
    # 入力形状を維持したまま結果を返す
    return eitp.(desired_pos_dim0, desired_pos_dim1)
end

function fast_trilinear_interpolate(
    pos0::AbstractArray,
    pos1::AbstractArray,
    pos2::AbstractArray,
    coords0::AbstractVector,
    coords1::AbstractVector,
    coords2::AbstractVector,
    values::AbstractArray{T,3},
) where {T}
    # 3次元のグリッド補間オブジェクトを作成
    # Gridded(Linear()) が Trilinear に相当します
    itp = interpolate((coords0, coords1, coords2), values, Gridded(Linear()))

    # 範囲外を NaN に設定
    sitp = extrapolate(itp, NaN)

    # ブロードキャストで一括計算
    return sitp.(pos0, pos1, pos2)
end

function fast_trilinear_interpolate_rectilinear(
    pos0,
    pos1,
    pos2,
    coords0,
    coords1,
    coords2,
    values,
)
    # 等間隔（Regular Grid）用の補間
    # 1. 単位格子(1, 2, 3...)上の補間として定義
    itp = interpolate(values, BSpline(Linear()))

    # 2. 実際の物理座標 (coords) にスケールを合わせる
    sitp = scale(itp, coords0, coords1, coords2)

    # 3. 範囲外を NaN に設定
    eitp = extrapolate(sitp, NaN)

    return eitp.(pos0, pos1, pos2)
end
