from copy import deepcopy
import unittest
from paired_stats_v1 import (invented_rows, validate_rows, paired_summary,
                             resampling_plan, cluster_bootstrap, exact_discordance_p)


class PairedMechanicsTests(unittest.TestCase):
    def setUp(self): self.rows = invented_rows()

    def test_hand_calculated_positive_stratum(self):
        value = paired_summary(self.rows, synthetic=True)
        self.assertEqual((value["n_total"], value["n_unknown_reference"], value["n_stratum"], value["n_clusters"]), (8, 1, 4, 3))
        self.assertEqual((value["candidate_rate"], value["comparator_rate"], value["paired_difference"]), (.75, .5, .25))
        self.assertEqual((value["gains"], value["losses"]), (2, 1))

    def test_hand_calculated_negative_stratum(self):
        value = paired_summary(self.rows, reference=0, synthetic=True)
        self.assertEqual(value["n_stratum"], 3)
        self.assertAlmostEqual(value["candidate_rate"], 1/3)
        self.assertAlmostEqual(value["comparator_rate"], 1/3)
        self.assertEqual(value["paired_difference"], 0)

    def test_repeated_cluster_draws_and_type7_endpoints(self):
        plan = [["c1", "c1", "c2"], ["c3", "c3", "c3"]]
        value = cluster_bootstrap(self.rows, plan, synthetic=True)
        self.assertEqual(value["values"], [.2, 1.0])
        self.assertAlmostEqual(value["lower"], .22)
        self.assertAlmostEqual(value["upper"], .98)

    def test_output_weighted_not_mean_of_cluster_means(self):
        plan = [["c1", "c2", "c3"], ["c1", "c2", "c3"]]
        result = cluster_bootstrap(self.rows, plan, synthetic=True)
        self.assertEqual(result["values"], [.25, .25])
        self.assertNotAlmostEqual(result["values"][0], (.5-1+1)/3)

    def test_shared_schedule_deterministic_and_input_order_invariant(self):
        first = resampling_plan(self.rows, draws=20, synthetic=True)
        self.assertEqual(first, resampling_plan(list(reversed(self.rows)), draws=20, synthetic=True))
        self.assertEqual(cluster_bootstrap(self.rows, first, synthetic=True), cluster_bootstrap(list(reversed(self.rows)), first, synthetic=True))

    def test_identical_method_outputs_have_zero_difference_and_interval(self):
        for row in self.rows: row["candidate"] = row["comparator"]
        self.assertEqual(paired_summary(self.rows, synthetic=True)["paired_difference"], 0)
        plan = resampling_plan(self.rows, draws=20, synthetic=True)
        result = cluster_bootstrap(self.rows, plan, synthetic=True)
        self.assertEqual((result["lower"], result["upper"]), (0., 0.))

    def test_exact_discordance_matches_independent_hand_example(self):
        rows = [dict(scope="SYNTHETIC", invented=True, item_id=str(i), cluster_id=str(i), reference=1, candidate=1, comparator=0) for i in range(4)]
        self.assertEqual(exact_discordance_p(rows, synthetic=True), .125)
        for row in rows: row["candidate"], row["comparator"] = row["comparator"], row["candidate"]
        self.assertEqual(exact_discordance_p(rows, synthetic=True), .125)
        for row in rows: row["candidate"] = row["comparator"]
        self.assertEqual(exact_discordance_p(rows, synthetic=True), 1.)

    def test_exact_test_rejects_cluster_dependence(self):
        with self.assertRaises(ValueError): exact_discordance_p(self.rows, synthetic=True)

    def test_bad_references_outputs_scopes_and_duplicates_rejected(self):
        for field, value in [("reference", True), ("reference", 2), ("candidate", None), ("candidate", True), ("comparator", -1), ("scope", "REAL"), ("invented", False)]:
            rows = deepcopy(self.rows); rows[0][field] = value
            with self.subTest(field=field, value=value), self.assertRaises(ValueError): validate_rows(rows, synthetic=True)
        rows = deepcopy(self.rows); rows[1]["item_id"] = rows[0]["item_id"]
        with self.assertRaises(ValueError): validate_rows(rows, synthetic=True)
        for marker in [False, 1]:
            with self.assertRaises(ValueError): paired_summary(self.rows, synthetic=marker)

    def test_absent_stratum_single_cluster_and_invalid_schedules_rejected(self):
        rows = deepcopy(self.rows)
        for row in rows: row["reference"] = None
        with self.assertRaises(ValueError): paired_summary(rows, synthetic=True)
        for row in self.rows: row["cluster_id"] = "only"
        with self.assertRaises(ValueError): cluster_bootstrap(self.rows, [["only"], ["only"]], synthetic=True)
        for plan in [[], [["c1"], ["c1"]], [["unknown"]*3]*2]:
            with self.assertRaises(ValueError): cluster_bootstrap(invented_rows(), plan, synthetic=True)

    def test_invented_fixture_never_mutated(self):
        original = deepcopy(self.rows)
        paired_summary(self.rows, synthetic=True)
        cluster_bootstrap(self.rows, resampling_plan(self.rows, draws=10, synthetic=True), synthetic=True)
        self.assertEqual(self.rows, original)


if __name__ == "__main__": unittest.main()
