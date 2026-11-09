"""One-dimensional bicycle model; Python 3.10+, standard library only.

Power is crank mechanical power; grade is rise/run; wind opposes along road.
Read 模型参数.json and 初稿.md for scenario definitions and parameter provenance.
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass, fields, replace
from itertools import product
import json
import math
from pathlib import Path
from typing import Iterable


def nonnegative(*values: float) -> None:
    if any(isinstance(x, bool) or not isinstance(x, (int, float))
           or not math.isfinite(x) or x < 0 for x in values):
        raise ValueError('Expected finite nonnegative numeric values.')


def positive(*values: float) -> None:
    nonnegative(*values)
    if any(x == 0 for x in values):
        raise ValueError('Expected strictly positive values.')


def teeth(*values: int) -> None:
    if any(type(x) is not int or x <= 0 for x in values):
        raise ValueError('Tooth counts must be positive integers.')


@dataclass(frozen=True)
class Model:
    total_mass_kg: float = 75.0
    circumference_m: float = 2.146
    rho: float = 1.225
    cda_m2: float = 0.30
    crr: float = 0.004
    eta: float = 0.97
    g: float = 9.80665
    rotational_equivalent_kg: float = 1.5

    def __post_init__(self) -> None:
        nonnegative(*self.__dict__.values())
        positive(self.total_mass_kg, self.circumference_m, self.rho, self.g, self.eta)
        if self.eta > 1:
            raise ValueError('Efficiency must be in (0, 1].')

    @property
    def mass(self) -> float:
        return self.total_mass_kg

    @property
    def effective_mass(self) -> float:
        return self.mass + self.rotational_equivalent_kg

    def forces(self, speed_kmh: float, grade: float = 0.0,
               headwind_kmh: float = 0.0) -> tuple[float, float, float]:
        nonnegative(speed_kmh, grade, headwind_kmh)
        v, w = speed_kmh / 3.6, headwind_kmh / 3.6
        denominator = math.hypot(1.0, grade)
        return (self.mass * self.g * grade / denominator,
                self.crr * self.mass * self.g / denominator,
                0.5 * self.rho * self.cda_m2 * (v + w) ** 2)

    def cruise(self, speed_kmh: float, grade: float = 0.0,
               headwind_kmh: float = 0.0) -> float:
        return sum(self.forces(speed_kmh, grade, headwind_kmh)) * speed_kmh / 3.6 / self.eta

    def speed_at_power(self, power_w: float, grade: float = 0.0,
                       headwind_kmh: float = 0.0) -> float:
        nonnegative(power_w, grade, headwind_kmh)
        if power_w == 0:
            return 0.0
        if self.cda_m2 == 0 and self.crr == 0 and grade == 0:
            raise ValueError('No finite positive-power steady speed without resistance.')
        lo, hi = 0.0, 20.0
        while self.cruise(hi, grade, headwind_kmh) < power_w:
            hi *= 2
            if hi > 1e6:
                raise ValueError('Steady speed is outside the supported search range.')
        for _ in range(90):
            mid = (lo + hi) / 2
            if self.cruise(mid, grade, headwind_kmh) < power_w:
                lo = mid
            else:
                hi = mid
        return (lo + hi) / 2

    def speed_at_cadence(self, front: int, rear: int, rpm: float) -> float:
        teeth(front, rear)
        nonnegative(rpm)
        return 0.06 * self.circumference_m * rpm * front / rear

    def cadence(self, front: int, rear: int, speed_kmh: float) -> float:
        teeth(front, rear)
        nonnegative(speed_kmh)
        return speed_kmh / (0.06 * self.circumference_m * front / rear)

    @staticmethod
    def crank_torque(power_w: float, rpm: float) -> float:
        nonnegative(power_w)
        positive(rpm)
        return power_w * 60 / (2 * math.pi * rpm)

    def constant_accel(self, start_kmh: float, end_kmh: float, seconds: float,
                       grade: float = 0.0, headwind_kmh: float = 0.0) -> dict[str, float]:
        nonnegative(start_kmh, end_kmh, grade, headwind_kmh)
        positive(seconds)
        if end_kmh <= start_kmh:
            raise ValueError('End speed must exceed start speed.')
        v0, v1, w = start_kmh / 3.6, end_kmh / 3.6, headwind_kmh / 3.6
        a = (v1 - v0) / seconds
        mean_v = (v0 + v1) / 2
        mean_v2 = (v0 * v0 + v0 * v1 + v1 * v1) / 3
        mean_v3 = (v0**3 + v0 * v0 * v1 + v0 * v1 * v1 + v1**3) / 4
        resistance = self.mass * self.g * (grade + self.crr) / math.hypot(1, grade)
        aero = 0.5 * self.rho * self.cda_m2
        mean_resistance = (resistance * mean_v + aero * (
            mean_v3 + 2 * w * mean_v2 + w * w * mean_v)) / self.eta
        kinetic = 0.5 * self.effective_mass * (v1 * v1 - v0 * v0)
        extra = kinetic / (self.eta * seconds)
        return {
            'acceleration_m_s2': a,
            'distance_m': seconds * mean_v,
            'kinetic_energy_increase_j': kinetic,
            'mean_resistance_power_w': mean_resistance,
            'mean_extra_accel_power_w': extra,
            'mean_total_power_w': mean_resistance + extra,
            'end_power_w': self.cruise(end_kmh, grade, headwind_kmh)
            + self.effective_mass * a * v1 / self.eta,
            'final_cruise_power_w': self.cruise(end_kmh, grade, headwind_kmh),
        }

    def constant_power_time(self, start_kmh: float, end_kmh: float, power_w: float,
                            grade: float = 0.0, headwind_kmh: float = 0.0) -> float:
        """Integrate dt/dv; reject unreachable/asymptotic targets and zero starts."""
        positive(start_kmh, end_kmh, power_w)
        nonnegative(grade, headwind_kmh)
        if end_kmh <= start_kmh:
            raise ValueError('End speed must exceed start speed.')
        if power_w <= self.cruise(end_kmh, grade, headwind_kmh):
            raise ValueError('Target is unreachable in finite time at this constant power.')
        lo, hi = start_kmh / 3.6, end_kmh / 3.6

        def integrand(v: float) -> float:
            return self.effective_mass * v / (self.eta * (
                power_w - self.cruise(v * 3.6, grade, headwind_kmh)))

        previous = None
        for exponent in range(6, 19):
            segments = 2**exponent
            step = (hi - lo) / segments
            value = (integrand(lo) + integrand(hi)
                     + 4 * sum(integrand(lo + i * step) for i in range(1, segments, 2))
                     + 2 * sum(integrand(lo + i * step) for i in range(2, segments, 2))) * step / 3
            if previous is not None and abs(value - previous) <= 1e-7 * max(1, abs(value)):
                return value
            previous = value
        raise ValueError('Time integral did not converge; target is too close to steady speed.')


def load_config(path: Path) -> tuple[Model, dict]:
    config = json.loads(path.read_text(encoding='utf-8-sig'))
    if not isinstance(config, dict) or set(config) != {'model', 'grid', 'aero_profiles', 'rolling_profiles'}:
        raise ValueError('Expected model, grid, aero_profiles and rolling_profiles sections.')
    if set(config['model']) != {field.name for field in fields(Model)}:
        raise ValueError('Model fields must match the documented schema; use total_mass_kg.')
    model = Model(**config['model'])
    grid = config['grid']
    required = {'mass_min_kg', 'mass_max_kg', 'mass_step_kg', 'grades_percent',
                'headwinds_kmh', 'speeds_kmh', 'powers_w', 'front_teeth',
                'cassette_teeth', 'cadence_min_rpm', 'cadence_max_rpm'}
    if set(grid) != required:
        raise ValueError('Unexpected or missing grid fields.')
    for key in ('mass_min_kg', 'mass_max_kg', 'mass_step_kg'):
        teeth(grid[key])
    if grid['mass_min_kg'] > grid['mass_max_kg'] or (
            grid['mass_max_kg'] - grid['mass_min_kg']) % grid['mass_step_kg']:
        raise ValueError('Mass endpoints must align with the mass step.')
    for key in ('grades_percent', 'headwinds_kmh', 'speeds_kmh', 'powers_w'):
        if not isinstance(grid[key], list) or not grid[key] or len(set(grid[key])) != len(grid[key]):
            raise ValueError('Grid lists must be nonempty and unique.')
        nonnegative(*grid[key])
    for key in ('speeds_kmh', 'powers_w'):
        positive(*grid[key])
    for key in ('front_teeth', 'cassette_teeth'):
        if not isinstance(grid[key], list) or not grid[key] or len(set(grid[key])) != len(grid[key]):
            raise ValueError('Gear lists must be nonempty and unique.')
        teeth(*grid[key])
    if grid['cassette_teeth'] != sorted(grid['cassette_teeth']):
        raise ValueError('Cassette teeth must be ordered from small to large.')
    positive(grid['cadence_min_rpm'], grid['cadence_max_rpm'])
    if grid['cadence_min_rpm'] >= grid['cadence_max_rpm']:
        raise ValueError('Cadence band must have increasing endpoints.')
    for section in ('aero_profiles', 'rolling_profiles'):
        profiles = config[section]
        if not isinstance(profiles, list) or not profiles:
            raise ValueError('Profiles must be a nonempty list.')
        if any(set(profile) != {'name', 'value', 'basis', 'source'} for profile in profiles):
            raise ValueError('Each profile needs name, value, basis and source.')
        positive(*(profile['value'] for profile in profiles))
        if len({profile['value'] for profile in profiles}) != len(profiles):
            raise ValueError('Profile values must be unique.')
        if any(not isinstance(profile[key], str) or not profile[key].strip()
               for profile in profiles for key in ('name', 'basis', 'source')):
            raise ValueError('Profile provenance must be nonempty text.')
    return model, config


def mass_levels(config: dict) -> range:
    grid = config['grid']
    return range(grid['mass_min_kg'], grid['mass_max_kg'] + 1, grid['mass_step_kg'])


def scenarios(model: Model, config: dict) -> Iterable[Model]:
    for mass, aero, rolling in product(mass_levels(config), config['aero_profiles'], config['rolling_profiles']):
        yield replace(model, total_mass_kg=mass, cda_m2=aero['value'], crr=rolling['value'])


def write_csv(path: Path, header: list[str], rows: Iterable) -> None:
    with path.open('w', encoding='utf-8-sig', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(header)
        for row in rows:
            writer.writerow([round(x, 6) if isinstance(x, float) else x for x in row])


def export_data(model: Model, config: dict, out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    grid = config['grid']

    def cruise_rows(models: Iterable[Model]) -> Iterable:
        for item in models:
            for grade, wind, speed in product(grid['grades_percent'], grid['headwinds_kmh'], grid['speeds_kmh']):
                power = item.cruise(speed, grade / 100, wind)
                yield [item.mass, item.cda_m2, item.crr, speed, grade, wind, power, power / item.mass]

    header = ['总质量_kg', 'CdA_m2', 'Crr', '地速_kmh', '坡度_percent', '正逆风_kmh', '曲柄功率_W', 'W每kg总质量']
    write_csv(out / '巡航功率完整数据.csv', header, cruise_rows([model]))
    write_csv(out / '质量分级巡航数据.csv', header, cruise_rows(scenarios(model, config)))
    sample = model.constant_accel(20, 35, 10)
    write_csv(out / '加速情景数据.csv', ['总质量_kg', 'CdA_m2', 'Crr', '起始速度_kmh', '目标速度_kmh', '加速时间_s', *sample],
              ([model.mass, model.cda_m2, model.crr, 20, 35, t, *model.constant_accel(20, 35, t).values()]
               for t in (10, 15, 30)))

    def climb_rows() -> Iterable:
        rear = max(grid['cassette_teeth'])
        for item in scenarios(model, config):
            for grade, power in product((6, 8, 10), grid['powers_w']):
                speed = item.speed_at_power(power, grade / 100)
                for front in grid['front_teeth']:
                    cadence = item.cadence(front, rear, speed)
                    yield [item.mass, item.cda_m2, item.crr, grade, 0, power, front, rear,
                           speed, cadence, item.crank_torque(power, cadence)]

    write_csv(out / '质量分级爬坡齿比数据.csv',
              ['总质量_kg', 'CdA_m2', 'Crr', '坡度_percent', '正逆风_kmh', '曲柄功率_W', '前盘_T', '后齿片_T', '稳态速度_kmh', '踏频_rpm', '曲柄扭矩_Nm'], climb_rows())

    def gear_rows() -> Iterable:
        for front, speed in product(grid['front_teeth'], (25, 30, 35, 40, 45)):
            for index, rear in enumerate(grid['cassette_teeth']):
                cadence = model.cadence(front, rear, speed)
                next_rear = grid['cassette_teeth'][index + 1] if index + 1 < len(grid['cassette_teeth']) else None
                yield [front, rear, speed, cadence,
                       grid['cadence_min_rpm'] <= cadence <= grid['cadence_max_rpm'],
                       100 * (next_rear / rear - 1) if next_rear else '',
                       model.speed_at_cadence(front, rear, grid['cadence_min_rpm']),
                       model.speed_at_cadence(front, rear, grid['cadence_max_rpm'])]

    write_csv(out / '常用速度档位数据.csv', ['前盘_T', '后齿片_T', '地速_kmh', '踏频_rpm', '位于设定踏频区间', '换到下一大齿片踏频增幅_percent', '最低踏频对应速度_kmh', '最高踏频对应速度_kmh'], gear_rows())

    def constant_power_rows() -> Iterable:
        for power in (150, 170, 201, 250, 300):
            reachable = power > model.cruise(35)
            duration = model.constant_power_time(20, 35, power) if reachable else ''
            yield [model.mass, model.cda_m2, model.crr, 20, 35, power,
                   model.speed_at_power(power), duration, '可达' if reachable else '不能有限时间达到']

    write_csv(out / '恒功率加速数据.csv', ['总质量_kg', 'CdA_m2', 'Crr', '起始速度_kmh', '目标速度_kmh', '曲柄功率_W', '稳态速度_kmh', '到达用时_s', '状态'], constant_power_rows())
    (out / '计算参数快照.json').write_text(json.dumps(config, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')


def run_checks(model: Model) -> None:
    """Checks never write data and remain active under python -O."""
    for grade, wind, power in product((0, .06, .10), (0, 10, 20), (150, 170, 210)):
        speed = model.speed_at_power(power, grade, wind)
        if abs(model.cruise(speed, grade, wind) - power) > 1e-8:
            raise ArithmeticError('Power inversion check failed.')
    steps, duration = 20000, 15
    acceleration = (35 - 20) / 3.6 / duration
    values = [model.cruise(20 + 15 * i / steps, .06, 10)
              + model.effective_mass * acceleration * ((20 + 15 * i / steps) / 3.6) / model.eta
              for i in range(steps + 1)]
    numeric = (values[0] / 2 + sum(values[1:-1]) + values[-1] / 2) / steps
    if abs(numeric - model.constant_accel(20, 35, duration, .06, 10)['mean_total_power_w']) > 1e-5:
        raise ArithmeticError('Acceleration integration check failed.')


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, default=Path(__file__).with_name('模型参数.json'))
    parser.add_argument('--output-dir', type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument('--check', action='store_true', help='Validate configuration and numerical checks without writing files.')
    args = parser.parse_args()
    model, config = load_config(args.config)
    run_checks(model)
    print('All numerical checks passed.')
    if not args.check:
        export_data(model, config, args.output_dir)
        print('Data saved in:', args.output_dir.resolve())


if __name__ == '__main__':
    main()
