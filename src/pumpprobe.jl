delaytime_fs(mirror_move_μm::Real) = 3.335640951981521 * mirror_move_μm

position_mm_to_delaytime_fs(position_mm::Real) = delaytime_fs(2 * position_mm * 1e3)
