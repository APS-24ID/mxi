"""mxeq truth: scaled observations against a refined structure's Fc^2.

On planted data -- P 21 21 21, a known scale, honest sigmas -- with a reflection
all of whose observations are 20 per cent low (as extinction would make them:
the model's part) and single observations 20 sigma off (as a processing fault
would: the processing's part): the scale comes back, the sigmas calibrate, the
single faults are the worst and their reflections' means stay near zero, and
the low reflection's mean is far below. And a .fcf in SHELXL's layout is read,
and a permuted cell's setting found.
"""

import numpy as np
import pytest

from mxeq import truth

FCF = """data_test
_cell_length_a 5.4220
_cell_length_b 8.1310
_cell_length_c 12.0230
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
loop_
 _space_group_symop_operation_xyz
 'x, y, z'
 '-x+1/2, -y, z+1/2'
 '-x, y+1/2, -z+1/2'
 'x+1/2, -y+1/2, -z'
loop_
 _refln_index_h
 _refln_index_k
 _refln_index_l
 _refln_F_squared_calc
 _refln_F_squared_meas
 _refln_F_squared_sigma
 _refln_observed_status
 1 1 1 120.5 118.0 2.1 o
 2 0 1 33.0 34.2 1.0 o
 0 2 3 7.5 7.1 0.6 o
"""


def test_an_fcf_is_read(tmp_path):
    path = tmp_path / "s.fcf"
    path.write_text(FCF)
    f = truth.read_fcf(str(path))
    assert f.cell[:3] == pytest.approx((5.422, 8.131, 12.023))
    assert f.name.replace(" ", "") == "P212121"
    assert f.hkl.tolist() == [[1, 1, 1], [2, 0, 1], [0, 2, 3]]
    assert f.fc2.tolist() == [120.5, 33.0, 7.5]
    assert f.sigma.tolist() == [2.1, 1.0, 0.6]


def test_a_permuted_cells_setting_is_found():
    op = truth.setting_from_cells(
        (12.023, 5.422, 8.131, 90, 90, 90), (5.422, 8.131, 12.023, 90, 90, 90)
    )
    assert op is not None and round(np.linalg.det(op)) == 1
    lengths = np.abs(op) @ np.array([12.023, 5.422, 8.131])
    assert lengths == pytest.approx([5.422, 8.131, 12.023])
    assert truth.parse_operator("b,c,a").tolist() == [[0, 1, 0], [0, 0, 1], [1, 0, 0]]


def planted(seed=3):
    rng = np.random.default_rng(seed)
    cell = (5.422, 8.131, 12.023, 90.0, 90.0, 90.0)
    hall = " P 2ac 2ab"  # P 21 21 21
    import gemmi

    ops = [np.array(op.rot) // op.DEN for op in gemmi.symops_from_hall(hall).sym_ops]
    unique = [(h, k, l) for h in range(0, 5) for k in range(0, 7) for l in range(1, 10)]
    fc2 = {u: float(rng.exponential(100.0)) + 1.0 for u in unique}
    hkl, i_obs, var, true_ref = [], [], [], []
    for u in unique:
        for r in ops:  # each reflection observed at its equivalents
            h = np.array(u) @ r
            sigma = np.sqrt(fc2[u] * 2.5 + 4.0)
            hkl.append(h)
            i_obs.append(2.5 * fc2[u] + rng.normal(0.0, sigma))
            var.append(sigma**2)
            true_ref.append(u)
    hkl, i_obs, var = np.array(hkl), np.array(i_obs), np.array(var)
    # The model's part: the strongest reflection, all its observations 20 per
    # cent low -- as extinction lowers the strong ones.
    low = max(unique, key=lambda u: fc2[u])
    for j, u in enumerate(true_ref):
        if u == low:
            i_obs[j] *= 0.8
    # The processing's part: single observations 20 sigma high.
    faults = [5, 41, 77]
    for j in faults:
        i_obs[j] += 20.0 * np.sqrt(var[j])
    fcf = truth.Fcf(
        cell=cell,
        hall=hall,
        hkl=np.array(unique),
        fc2=np.array([fc2[u] for u in unique]),
        fo2=np.zeros(len(unique)),
        sigma=np.ones(len(unique)),
        name="P 21 21 21",
    )
    obs = truth.Observations(
        hkl=hkl, intensity=i_obs, variance=var, cell=cell, hall=hall
    )
    return obs, fcf, low, faults, true_ref


def test_the_scale_comes_back_and_the_sigmas_calibrate():
    obs, fcf, *_ = planted()
    res = truth.analyse(obs, fcf, anomalous=False, bins=4)
    assert res.matched == len(obs.hkl) and res.unmatched == 0
    assert res.k == pytest.approx(2.5, rel=0.01)
    assert res.robust_var_z == pytest.approx(1.0, abs=0.15)


def test_a_single_fault_is_the_processings_and_a_shared_offset_the_models():
    obs, fcf, low, faults, true_ref = planted()
    res = truth.analyse(obs, fcf, anomalous=False, bins=4, worst=10)
    worst_rows = [int(res.rows[np.argsort(-np.abs(res.z))[i]]) for i in range(3)]
    assert sorted(worst_rows) == sorted(faults)
    for j in faults:  # the reflection's own mean stays far nearer zero
        assert abs(res.zbar[j]) < abs(res.z[j]) / 2.5
    lows = [j for j, u in enumerate(true_ref) if u == low]
    assert all(res.zbar[j] < -3.0 for j in lows)  # all its observations low together


def test_a_permuted_setting_matches_every_observation():
    obs, fcf, *_ = planted()
    # The observations indexed in the order c,a,b of the .fcf's axes.
    to_obs = truth.parse_operator("c,a,b")
    obs.hkl = obs.hkl @ to_obs.T
    obs.cell = (12.023, 5.422, 8.131, 90.0, 90.0, 90.0)
    op = truth.setting_from_cells(obs.cell, fcf.cell)
    res = truth.analyse(obs, fcf, reindex=op, anomalous=False, bins=4)
    assert res.unmatched == 0
    assert res.k == pytest.approx(2.5, rel=0.01)
