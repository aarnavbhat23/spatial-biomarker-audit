import itertools
import unittest
import numpy as np
from joint_measurement import enumerate_worlds, measurement_actions, pair_requirements, optimal_fixed_plan


def transcript_signs(worlds, actions, selected, error):
    records = {}
    for world in worlds:
        options = []
        for i in selected:
            action = actions[i]
            key = 'fragile_present' if action['kind'] == 'geometry' else 'labels'
            true = world[key][action['index']]
            options.append(range(3) if action['kind'] == 'routine' else (true,))
        for answers in itertools.product(*options):
            errors = sum(answer != world['labels'][actions[i]['index']]
                         for i, answer in zip(selected, answers) if actions[i]['kind'] == 'routine')
            if errors <= error:
                records.setdefault(answers, set()).add(world['sign'])
    return records


class JointMeasurementTests(unittest.TestCase):
    def test_pair_criterion_matches_every_transcript_for_every_plan(self):
        worlds = enumerate_worlds((1, 2, 0), ((0, 1), (1, 2)), ((0, 2),))
        actions = measurement_actions(3, ((0, 2),))
        for error in (0, 1):
            _, weights, quota = pair_requirements(worlds, actions, error)
            for mask in range(1 << len(actions)):
                selected = [i for i in range(len(actions)) if mask & (1 << i)]
                predicted = bool(np.all(weights[:, selected].sum(axis=1) >= quota))
                actual = all(len(signs) == 1 for signs in transcript_signs(worlds, actions, selected, error).values())
                self.assertEqual(predicted, actual, (error, mask))

    def test_exact_measurements_always_suffice_and_errors_cannot_help(self):
        worlds = enumerate_worlds((1, 2, 0), ((0, 1), (1, 2)), ((0, 2),))
        actions = measurement_actions(3, ((0, 2),))
        exact = optimal_fixed_plan(worlds, actions, 1, {'adjudication', 'geometry'})
        zero = optimal_fixed_plan(worlds, actions, 0)
        noisy = optimal_fixed_plan(worlds, actions, 1)
        self.assertIsNotNone(exact['cost'])
        self.assertLessEqual(zero['cost'], noisy['cost'])
        self.assertLessEqual(noisy['cost'], exact['cost'])

    def test_enumeration_respects_exclusivity_and_complete_product(self):
        worlds = enumerate_worlds((1, 2, 0), ((0, 1), (1, 2)), ((0, 2),), (3, 3))
        self.assertEqual(len(worlds), 3**3 * 2)
        self.assertEqual(len({(tuple(w['labels']), tuple(map(tuple, w['edges']))) for w in worlds}), len(worlds))
        for w in worlds:
            self.assertEqual(w['sign'], (w['numerator'] > 0) - (w['numerator'] < 0))

    def test_invalid_graph_and_error_budgets_fail(self):
        with self.assertRaises(ValueError):
            enumerate_worlds((1, 2), ((0, 1),), ((0, 1),))
        with self.assertRaises(ValueError):
            pair_requirements([], [], -.1)


if __name__ == '__main__':
    unittest.main()
