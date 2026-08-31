import unittest

import numpy as np

from operaattori.credit import (
    CreditConfig,
    evaluate_ports,
    selectivity,
    train,
)


class CreditTests(unittest.TestCase):
    def test_selectivity_bounds(self):
        self.assertAlmostEqual(selectivity(1.0, 0.0), 1.0)
        self.assertLess(selectivity(0.0, 1.0), -0.999999)

    def test_no_credit_does_not_grow_large_material(self):
        c = CreditConfig(trials=8)
        row = train("no_credit", 0, c)
        self.assertLess(float(np.sum(row["morphology"])), 2.0)

    def test_eval_is_deterministic_for_frozen_morphology(self):
        c = CreditConfig(eval_steps=120, eval_pulse_steps=20)
        m = np.linspace(0.0, 0.2, c.n)
        a0 = evaluate_ports(m, c)
        a1 = evaluate_ports(m, c)
        self.assertEqual(a0, a1)


if __name__ == "__main__":
    unittest.main()
