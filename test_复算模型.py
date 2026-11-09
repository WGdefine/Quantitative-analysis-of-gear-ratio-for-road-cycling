"""Numerical, domain and configuration regression checks (standard-library unittest)."""
import csv
from dataclasses import replace
import json
import math
from pathlib import Path
import tempfile
import unittest

from 复算模型 import Model, export_data, load_config, mass_levels


class ModelTests(unittest.TestCase):
    def setUp(self):
        self.model = Model()

    def test_original_reference_regression(self):
        old = Model(total_mass_kg=73, cda_m2=.32)
        self.assertAlmostEqual(old.cruise(35), 214.38779291997227, places=9)
        self.assertAlmostEqual(old.constant_accel(20, 35, 15)['mean_total_power_w'], 282.2907672899608, places=9)

    def test_power_matches_independent_angle_formula(self):
        m = self.model
        for mass in (50, 75, 120):
            item = replace(m, total_mass_kg=mass)
            for grade in (0, .06, .1):
                for wind in (0, 20):
                    v, w, theta = 35 / 3.6, wind / 3.6, math.atan(grade)
                    expected = (mass * m.g * (math.sin(theta) + m.crr * math.cos(theta))
                                + .5 * m.rho * m.cda_m2 * (v + w)**2) * v / m.eta
                    self.assertAlmostEqual(item.cruise(35, grade, wind), expected, places=10)

    def test_mass_interpolation_and_rotating_mass_scope(self):
        m = self.model
        a, b = replace(m, total_mass_kg=70), replace(m, total_mass_kg=75)
        exact = replace(m, total_mass_kg=73)
        self.assertAlmostEqual(exact.cruise(35, .06), .4 * a.cruise(35, .06) + .6 * b.cruise(35, .06), places=10)
        spin = replace(m, rotational_equivalent_kg=3)
        self.assertEqual(spin.cruise(35, .08), m.cruise(35, .08))
        self.assertGreater(spin.constant_accel(20, 35, 15)['mean_total_power_w'], m.constant_accel(20, 35, 15)['mean_total_power_w'])

    def test_pure_aero_cubic_and_wind_ground_speed(self):
        m = replace(self.model, crr=0)
        self.assertAlmostEqual(m.cruise(40), 8 * m.cruise(20), places=10)
        self.assertAlmostEqual(m.cruise(20, 0, 20), 4 * m.cruise(20), places=10)
        self.assertEqual(m.cruise(0, .1, 20), 0)

    def test_inverse_and_gearing_round_trip(self):
        for mass in (50, 75, 120):
            m = replace(self.model, total_mass_kg=mass)
            for power in (100, 170, 300):
                for grade in (0, .08):
                    speed = m.speed_at_power(power, grade, 10)
                    self.assertAlmostEqual(m.cruise(speed, grade, 10), power, places=9)
            self.assertAlmostEqual(m.cadence(44, 17, m.speed_at_cadence(44, 17, 90)), 90)

    def test_acceleration_energy_balance_independent_quadrature(self):
        m = self.model
        n, t, v0, v1 = 10000, 15, 20 / 3.6, 35 / 3.6
        theta, wind = math.atan(.06), 10 / 3.6
        def power(v):
            forces = m.mass * m.g * (math.sin(theta) + m.crr * math.cos(theta)) + .5 * m.rho * m.cda_m2 * (v + wind)**2
            return (forces + m.effective_mass * (v1 - v0) / t) * v / m.eta
        values = [power(v0 + (v1 - v0) * i / n) for i in range(n + 1)]
        numeric = (values[0] + values[-1] + 4 * sum(values[1:n:2]) + 2 * sum(values[2:n:2])) / (3 * n)
        self.assertAlmostEqual(numeric, m.constant_accel(20, 35, t, .06, 10)['mean_total_power_w'], places=8)

    def test_constant_power_time_with_independent_rk4(self):
        m, power = self.model, 250
        duration = m.constant_power_time(20, 35, power)
        dt, v = duration / 4000, 20 / 3.6
        def acceleration(v):
            resistance = m.mass * m.g * m.crr + .5 * m.rho * m.cda_m2 * v * v
            return (m.eta * power / v - resistance) / m.effective_mass
        for _ in range(4000):
            a = acceleration(v)
            b = acceleration(v + dt * a / 2)
            c = acceleration(v + dt * b / 2)
            d = acceleration(v + dt * c)
            v += dt * (a + 2 * b + 2 * c + d) / 6
        self.assertAlmostEqual(v * 3.6, 35, places=6)
        for power in (150, m.cruise(35)):
            with self.assertRaises(ValueError):
                m.constant_power_time(20, 35, power)
        with self.assertRaises(ValueError):
            m.constant_power_time(0, 35, 300)

    def test_invalid_inputs(self):
        m = self.model
        for invalid in (float('nan'), float('inf'), -1, True, '90'):
            for call in (lambda: m.speed_at_cadence(44, 11, invalid),
                         lambda: m.crank_torque(170, invalid),
                         lambda: m.cruise(invalid),
                         lambda: Model(total_mass_kg=invalid)):
                with self.assertRaises(ValueError):
                    call()
        for invalid in (0, -1, 44.5, True, float('nan')):
            with self.assertRaises(ValueError):
                m.cadence(invalid, 34, 20)
        with self.assertRaises(ValueError):
            m.crank_torque(170, 0)

    def test_configuration_controls_outputs_and_is_not_overwritten(self):
        path = Path(__file__).with_name('模型参数.json')
        _, config = load_config(path)
        self.assertEqual(list(mass_levels(config)), list(range(50, 121, 5)))
        config['model']['total_mass_kg'] = 80
        config['grid']['mass_min_kg'] = 70
        config['grid']['mass_max_kg'] = 75
        config['grid']['speeds_kmh'] = [20, 35]
        with tempfile.TemporaryDirectory() as directory:
            out = Path(directory)
            source = out / 'custom.json'
            source.write_text(json.dumps(config), encoding='utf-8')
            before = source.read_bytes()
            m, loaded = load_config(source)
            export_data(m, loaded, out)
            self.assertEqual(source.read_bytes(), before)
            with (out / '质量分级巡航数据.csv').open(encoding='utf-8-sig', newline='') as file:
                rows = list(csv.DictReader(file))
            self.assertEqual(len(rows), 2 * 3 * 3 * 6 * 3 * 2)
            self.assertEqual(len({tuple(row[k] for k in list(row)[:6]) for row in rows}), len(rows))
            self.assertNotIn('W每kg骑手', rows[0])
            with (out / '巡航功率完整数据.csv').open(encoding='utf-8-sig', newline='') as file:
                row = next(csv.DictReader(file))
            self.assertEqual(float(row['总质量_kg']), 80)


if __name__ == '__main__':
    unittest.main()
