"""Reproducible one-dimensional bicycle gearing / power model.

All values are scenario assumptions, not rider- or equipment-specific tests.
References for physical laws: see accompanying article.
Python 3.10+, standard library only. No network access or third-party dependencies.
Grade is a fraction: 6% = 0.06. Wind is a nonnegative opposing component ALONG road.
Power is crank mechanical power. Model excludes wheel aerodynamic spin torque,
separate hub-bearing drag, rider limb internal work, changing drivetrain losses.
"""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import csv
import json
import math

@dataclass(frozen=True)
class Model:
    rider_kg: float = 63.0
    bicycle_and_kit_kg: float = 10.0
    circumference_m: float = 2.146
    rho: float = 1.225
    cda_m2: float = 0.32
    crr: float = 0.004
    eta: float = 0.97
    g: float = 9.80665
    rotational_equivalent_kg: float = 1.5

    def __post_init__(self) -> None:
        values = self.__dict__.values()
        if not all(math.isfinite(x) for x in values):
            raise ValueError('All parameters must be finite.')
        if min(self.rider_kg, self.circumference_m, self.rho, self.g) <= 0:
            raise ValueError('Mass, circumference, density and gravity must be positive.')
        if min(self.bicycle_and_kit_kg, self.cda_m2, self.crr, self.rotational_equivalent_kg) < 0:
            raise ValueError('Mass additions and resistance parameters must be nonnegative.')
        if not 0 < self.eta <= 1:
            raise ValueError('Efficiency must be in (0, 1].')

    @property
    def mass(self) -> float:
        return self.rider_kg + self.bicycle_and_kit_kg

    @property
    def effective_mass(self) -> float:
        return self.mass + self.rotational_equivalent_kg

    @staticmethod
    def _check(*vals: float) -> None:
        if any(not math.isfinite(x) or x < 0 for x in vals):
            raise ValueError('This model is limited to finite nonnegative speed, grade and headwind.')

    def forces(self, speed_kmh: float, grade: float = 0.0, headwind_kmh: float = 0.0) -> tuple[float, float, float]:
        self._check(speed_kmh, grade, headwind_kmh)
        v, w = speed_kmh / 3.6, headwind_kmh / 3.6
        den = math.sqrt(1.0 + grade * grade)
        gravity = self.mass * self.g * grade / den
        rolling = self.crr * self.mass * self.g / den
        aero = 0.5 * self.rho * self.cda_m2 * (v + w) ** 2
        return gravity, rolling, aero

    def cruise(self, speed_kmh: float, grade: float = 0.0, headwind_kmh: float = 0.0) -> float:
        return sum(self.forces(speed_kmh, grade, headwind_kmh)) * speed_kmh / 3.6 / self.eta

    def speed_at_power(self, power_w: float, grade: float = 0.0, headwind_kmh: float = 0.0) -> float:
        self._check(power_w, grade, headwind_kmh)
        if power_w == 0:
            return 0.0
        if self.cda_m2 == 0 and self.crr == 0 and grade == 0:
            raise ValueError('No finite positive-power steady speed without any resistance.')
        lo, hi = 0.0, 20.0
        while self.cruise(hi, grade, headwind_kmh) < power_w:
            hi *= 2.0
            if hi > 1e6:
                raise ValueError('No suitable speed root in model range.')
        for _ in range(90):
            mid = (lo + hi) / 2.0
            if self.cruise(mid, grade, headwind_kmh) < power_w:
                lo = mid
            else:
                hi = mid
        return (lo + hi) / 2.0

    def speed_at_cadence(self, front: int, rear: int, rpm: float) -> float:
        if front <= 0 or rear <= 0 or rpm < 0:
            raise ValueError('Tooth counts must be positive and cadence nonnegative.')
        return 0.06 * self.circumference_m * rpm * front / rear

    def cadence(self, front: int, rear: int, speed_kmh: float) -> float:
        self._check(speed_kmh)
        if front <= 0 or rear <= 0:
            raise ValueError('Tooth counts must be positive.')
        return speed_kmh / (0.06 * self.circumference_m * front / rear)

    @staticmethod
    def crank_torque(power_w: float, rpm: float) -> float:
        if power_w < 0 or rpm <= 0:
            raise ValueError('Use nonnegative power and strictly positive cadence.')
        return power_w * 60.0 / (2.0 * math.pi * rpm)

    def constant_accel(self, start_kmh: float, end_kmh: float, seconds: float,
                       grade: float = 0.0, headwind_kmh: float = 0.0) -> dict[str, float]:
        self._check(start_kmh, end_kmh, grade, headwind_kmh)
        if not math.isfinite(seconds) or seconds <= 0 or end_kmh <= start_kmh:
            raise ValueError('Positive finite time and strictly increasing speed required.')
        v0, v1, w = start_kmh / 3.6, end_kmh / 3.6, headwind_kmh / 3.6
        a = (v1 - v0) / seconds
        mean_v = (v0 + v1) / 2.0
        mean_v2 = (v0*v0 + v0*v1 + v1*v1) / 3.0
        mean_v3 = (v0**3 + v0*v0*v1 + v0*v1*v1 + v1**3) / 4.0
        A = self.mass * self.g * (grade + self.crr) / math.sqrt(1 + grade*grade)
        B = 0.5 * self.rho * self.cda_m2
        mean_resistance_power = (A*mean_v + B*(mean_v3 + 2*w*mean_v2 + w*w*mean_v)) / self.eta
        kinetic_j = 0.5 * self.effective_mass * (v1*v1 - v0*v0)
        mean_accel_power = kinetic_j / (self.eta * seconds)
        return {
            'acceleration_m_s2': a,
            'distance_m': seconds * mean_v,
            'kinetic_energy_increase_j': kinetic_j,
            'mean_resistance_power_w': mean_resistance_power,
            'mean_extra_accel_power_w': mean_accel_power,
            'mean_total_power_w': mean_resistance_power + mean_accel_power,
            'end_power_w': self.cruise(end_kmh, grade, headwind_kmh) + self.effective_mass*a*v1/self.eta,
            'final_cruise_power_w': self.cruise(end_kmh, grade, headwind_kmh),
        }


def main() -> None:
    model = Model()
    out = Path(__file__).resolve().parent
    with (out / '巡航功率完整数据.csv').open('w', encoding='utf-8-sig', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['地速_kmh', '坡度_percent', '正逆风_kmh', '曲柄功率_W', 'W每kg总质量', 'W每kg骑手'])
        for grade_pct in (0, 2, 4, 6, 8, 10):
            for wind in (0, 10, 20):
                for speed in range(5, 51):
                    p = model.cruise(speed, grade_pct / 100, wind)
                    writer.writerow([speed, grade_pct, wind, round(p, 6), round(p/model.mass, 6), round(p/model.rider_kg, 6)])
    with (out / '加速情景数据.csv').open('w', encoding='utf-8-sig', newline='') as f:
        writer = csv.writer(f)
        sample = model.constant_accel(20, 35, 10)
        writer.writerow(['起始速度_kmh','目标速度_kmh','加速时间_s', *sample.keys()])
        for seconds in (10,15,30):
            result = model.constant_accel(20,35,seconds)
            writer.writerow([20,35,seconds,*[round(x,6) for x in result.values()]])
    (out / '模型参数.json').write_text(json.dumps(model.__dict__, ensure_ascii=False, indent=2),encoding='utf-8')
    # Independent checks: invert the steady model; zero-ramp limit; numerical integration of variable power.
    for grade in (0.0,0.06,0.10):
        for wind in (0.0,10.0,20.0):
            v = model.speed_at_power(170,grade,wind)
            assert abs(model.cruise(v,grade,wind)-170)<1e-8
    N=20000
    T=15
    a=(35-20)/3.6/T
    power=[]
    for i in range(N+1):
        speed=20+(35-20)*i/N
        power.append(model.cruise(speed,0.06,10)+model.effective_mass*a*(speed/3.6)/model.eta)
    numeric=(0.5*power[0]+sum(power[1:-1])+0.5*power[-1])/N
    exact=model.constant_accel(20,35,T,0.06,10)['mean_total_power_w']
    assert abs(numeric-exact)<1e-5
    print('All numerical checks passed. Data saved in:', out)
    print('\nCruise grades: rows 0,2,4,6,8,10%; columns 10,15,20,25,30,35 km/h')
    for pct in (0,2,4,6,8,10):
        print(pct,[round(model.cruise(v,pct/100)) for v in (10,15,20,25,30,35)])
    print('\nWind flat: columns 0,10,20 km/h; rows 20..45')
    for speed in (20,25,30,35,40,45):
        print(speed,[round(model.cruise(speed,0,w)) for w in (0,10,20)])
    print('\nAccel 20 -> 35')
    for t in (10,15,30): print(t,model.constant_accel(20,35,t))
    print('\n170 W inverse at 6% and 8%')
    for grade in (.06,.08):
        speed=model.speed_at_power(170,grade)
        print(grade, speed,[(f,round(model.cadence(f,34,speed),2),round(model.crank_torque(170,model.cadence(f,34,speed)),2)) for f in (50,46,44,42)])
    print('\n35 flat split:',[x*(35/3.6)/model.eta for x in model.forces(35)])
    print('Sensitivity35 CdA delta .04:',.5*model.rho*.04*(35/3.6)**3/model.eta)
    print('Sensitivity35 Crr delta .001:',model.mass*model.g*.001*(35/3.6)/model.eta)
    print('Cadence@35:50x16,46x15,44x14', [model.cadence(f,r,35) for f,r in ((50,16),(46,15),(44,14))])
    print('Gear@90 speeds',[(f,model.speed_at_cadence(f,11,90)) for f in (50,48,46,44,42)])

if __name__=='__main__':
    main()
