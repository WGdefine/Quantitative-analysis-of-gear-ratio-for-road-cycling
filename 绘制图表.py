"""Reproduce publication figures from the same JSON as the numerical model.

Optional dependency: matplotlib (only this plotting script needs it).
"""
from __future__ import annotations

import argparse
from dataclasses import replace
import os
from pathlib import Path
import tempfile

from 复算模型 import Model, load_config, mass_levels, run_checks


def draw_figures(model: Model, config: dict, out: Path) -> None:
    os.environ.setdefault('MPLCONFIGDIR', str(Path(tempfile.gettempdir()) / 'cycling_matplotlib'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap

    plt.rcParams.update({
        'font.sans-serif': ['Microsoft YaHei', 'SimHei', 'Noto Sans CJK SC', 'DejaVu Sans'],
        'axes.unicode_minus': False, 'font.size': 11,
        'axes.spines.top': False, 'axes.spines.right': False,
        'axes.titlepad': 16, 'savefig.facecolor': 'white',
    })
    out.mkdir(parents=True, exist_ok=True)
    colors = ['#2266A5', '#D17827', '#26856C', '#9D4E6A', '#7864A4', '#6E675D']
    note = (f'模型计算；总质量{model.mass:g}kg；CdA={model.cda_m2:g}m²；'
            f'Crr={model.crr:g}；ρ={model.rho:g}kg/m³；η={model.eta:g}')

    def save(fig, name: str, footer: str = note) -> None:
        fig.text(.5, .025, footer, ha='center', fontsize=10, color='#444444')
        fig.tight_layout(rect=(0, .065, 1, 1))
        fig.savefig(out / name, dpi=160)
        plt.close(fig)

    speeds = [5 + i / 5 for i in range(176)]
    fig, ax = plt.subplots(figsize=(12, 7))
    for index, grade in enumerate(config['grid']['grades_percent']):
        ax.plot(speeds, [model.cruise(v, grade / 100) for v in speeds],
                color=colors[index % len(colors)], label=f'{grade:g}%坡', linewidth=2)
    ax.set(title='图1  无风巡航：坡度如何改变功率需求',
           xlabel='地面速度（km/h）', ylabel='曲柄机械功率（W）', ylim=(0, None))
    ax.grid(alpha=.2)
    ax.legend(loc='upper left')
    save(fig, '图1_匀速巡航_坡度.png')

    fig, ax = plt.subplots(figsize=(12, 7))
    for index, wind in enumerate(config['grid']['headwinds_kmh']):
        ax.plot(speeds, [model.cruise(v, 0, wind) for v in speeds],
                color=colors[index % len(colors)], label=f'正面逆风{wind:g}km/h', linewidth=2)
    ax.set(title='图2  平路巡航：正面逆风如何改变功率需求',
           xlabel='地面速度（km/h）', ylabel='曲柄机械功率（W）', ylim=(0, None))
    ax.grid(alpha=.2)
    ax.legend(loc='upper left')
    save(fig, '图2_匀速巡航_逆风.png')

    fig, ax = plt.subplots(figsize=(12, 7))
    final_power = model.cruise(35)
    ax.axhline(final_power, color='#444444', linestyle='--', linewidth=1.5,
               label=f'35km/h巡航：{final_power:.1f}W')
    for index, duration in enumerate((10, 15, 30)):
        result = model.constant_accel(20, 35, duration)
        times = [duration * i / 200 for i in range(201)]
        acceleration = (35 - 20) / 3.6 / duration
        powers = [model.cruise(20 + 15 * t / duration)
                  + model.effective_mass * acceleration * ((20 + 15 * t / duration) / 3.6) / model.eta
                  for t in times]
        ax.plot(times, powers, color=colors[index], linewidth=2.2,
                label=f'{duration}秒匀加速（过程平均{result["mean_total_power_w"]:.0f}W）')
        ax.plot([duration, duration, 34], [powers[-1], final_power, final_power],
                color=colors[index], linestyle=':', linewidth=2)
        ax.scatter([duration], [powers[-1]], color=colors[index], zorder=4)
        ax.annotate(f'{powers[-1]:.0f}W', (duration, powers[-1]), xytext=(6, 8),
                    textcoords='offset points', color=colors[index])
    ax.set(title='图3  从20到35km/h：匀加速结束后转入巡航',
           xlabel='加速开始后的时间（s）', ylabel='瞬时曲柄机械功率（W）', xlim=(0, 35), ylim=(0, None))
    ax.grid(alpha=.2)
    ax.legend(loc='upper right', fontsize=10)
    save(fig, '图3_加速与巡航.png', note + '\n无风平路；虚线段表示加速结束；不模拟人体动作的平滑过渡')

    masses = list(mass_levels(config))
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    for index, profile in enumerate(config['aero_profiles']):
        models = [replace(model, total_mass_kg=mass, cda_m2=profile['value']) for mass in masses]
        label = f'CdA={profile["value"]:g}m²'
        axes[0].plot(masses, [item.cruise(35) for item in models], color=colors[index], label=label, marker='o', markersize=3)
        axes[1].plot(masses, [item.speed_at_power(170, .08) for item in models], color=colors[index], label=label, marker='o', markersize=3)
    axes[0].set(title='平路35km/h所需功率', xlabel='人车装备总质量（kg）', ylabel='曲柄机械功率（W）')
    axes[1].set(title='8%坡、170W可维持速度', xlabel='人车装备总质量（kg）', ylabel='稳态地速（km/h）')
    for ax in axes:
        ax.grid(alpha=.2)
        ax.legend()
    fig.suptitle('图4  总质量与风阻分别改变什么？', fontsize=16)
    save(fig, '图4_质量与风阻分级.png', f'模型计算；无风；Crr={model.crr:g}；其他参数见配置。质量与CdA独立取值，不按体重自动匹配风阻。')

    fronts = config['grid']['front_teeth']
    rears = config['grid']['cassette_teeth']
    low, high = config['grid']['cadence_min_rpm'], config['grid']['cadence_max_rpm']
    rpms = [[model.cadence(front, rear, 35) for rear in rears] for front in fronts]
    flags = [[int(low <= rpm <= high) for rpm in row] for row in rpms]
    fig, ax = plt.subplots(figsize=(13, 5.5))
    ax.imshow(flags, cmap=ListedColormap(['#F0F1F2', '#B4DDC9']), vmin=0, vmax=1, aspect='auto')
    for row, values in enumerate(rpms):
        for col, rpm in enumerate(values):
            ax.text(col, row, f'{rpm:.1f}', ha='center', va='center', fontsize=11,
                    color='#174A35' if flags[row][col] else '#777777')
    ax.set_xticks(range(len(rears)), [str(rear) for rear in rears])
    ax.set_yticks(range(len(fronts)), [str(front) for front in fronts])
    ax.set(title=f'图5  35km/h各档踏频：绿色为{low:g}–{high:g}rpm',
           xlabel='后齿片齿数（T）', ylabel='前盘齿数（T）')
    save(fig, '图5_常用巡航档位.png', f'纯传动几何计算；轮周{model.circumference_m:g}m；格内单位rpm；飞轮齿数组合为配置示例，装车兼容性需另核实。')


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, default=Path(__file__).with_name('模型参数.json'))
    parser.add_argument('--output-dir', type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    model, config = load_config(args.config)
    run_checks(model)
    draw_figures(model, config, args.output_dir)
    print('Five figures saved in:', args.output_dir.resolve())


if __name__ == '__main__':
    main()
