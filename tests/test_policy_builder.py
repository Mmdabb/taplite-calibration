from __future__ import annotations

import numpy as np

from taplite_calibration.policy_builder import _path_screen_pairs


def test_path_screen_incidence_counts_each_member_link() -> None:
    # Path 0 traverses two links in screen 10 and one link in screen 20.
    # Path 1 traverses one link outside the selected screens.
    links = np.array([1, 2, 3, 4], dtype=np.int32)
    link_offsets = np.array([0, 3, 4], dtype=np.int32)
    screen_lookup = np.array([-1, 0, 0, 1, -1], dtype=np.int32)

    path, screen, multiplicity = _path_screen_pairs(
        links, link_offsets, screen_lookup, screen_count=2
    )

    assert path.tolist() == [0, 0]
    assert screen.tolist() == [0, 1]
    assert multiplicity.tolist() == [2, 1]
