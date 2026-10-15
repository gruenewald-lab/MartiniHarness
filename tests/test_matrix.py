from martini_harness.matrix import expand_matrix


def test_product_of_axes():
    combos = expand_matrix({"target": ["a", "b"], "solvent": ["W", "CLF"]})
    assert combos == [
        {"target": "a", "solvent": "W"},
        {"target": "a", "solvent": "CLF"},
        {"target": "b", "solvent": "W"},
        {"target": "b", "solvent": "CLF"},
    ]


def test_exclude_partial_match():
    combos = expand_matrix(
        {"target": ["a", "b"], "solvent": ["W", "CLF"]},
        exclude=[{"solvent": "CLF", "target": "b"}, {"target": "a", "solvent": "W"}],
    )
    assert combos == [{"target": "a", "solvent": "CLF"}, {"target": "b", "solvent": "W"}]


def test_include_extends_matching_combinations():
    combos = expand_matrix(
        {"target": ["a", "b"], "solvent": ["W", "CLF"]},
        include=[{"solvent": "W", "density": 997}],
    )
    assert {"target": "a", "solvent": "W", "density": 997} in combos
    assert {"target": "b", "solvent": "W", "density": 997} in combos
    assert {"target": "a", "solvent": "CLF"} in combos
    assert len(combos) == 4


def test_include_without_axis_keys_extends_all():
    combos = expand_matrix({"solvent": ["W", "CLF"]}, include=[{"n_solvent": 10}])
    assert all(combo["n_solvent"] == 10 for combo in combos)


def test_include_never_overwrites_original_values_but_adds_new_combo():
    combos = expand_matrix({"solvent": ["W"]}, include=[{"solvent": "OCO", "density": 1}])
    assert combos == [{"solvent": "W"}, {"solvent": "OCO", "density": 1}]


def test_include_may_overwrite_added_values():
    combos = expand_matrix(
        {"solvent": ["W", "CLF"]},
        include=[{"density": 900}, {"solvent": "W", "density": 997}],
    )
    assert combos == [{"solvent": "W", "density": 997}, {"solvent": "CLF", "density": 900}]


def test_include_only_matrix():
    combos = expand_matrix({}, include=[{"target": "a", "solvent": "W"}, {"target": "b", "solvent": "CLF"}])
    assert combos == [{"target": "a", "solvent": "W"}, {"target": "b", "solvent": "CLF"}]


def test_exclude_happens_before_include():
    combos = expand_matrix(
        {"solvent": ["W", "CLF"]},
        exclude=[{"solvent": "W"}],
        include=[{"solvent": "W", "density": 997}],
    )
    assert combos == [{"solvent": "CLF"}, {"solvent": "W", "density": 997}]
