# Repository Guidelines

## Project Structure & Module Organization

This repository is a reproducible road-bicycle gearing and power study. All files live at the repository root:

- `初稿.md`: main Chinese article, numbered equations, tables, references, and relative image links.
- `计算附录.md`: acceleration derivations and numerical reference values.
- `复算模型.py`: Python model, data generation, and numerical checks.
- `模型参数.json`: exported default parameters.
- `巡航功率完整数据.csv` and `加速情景数据.csv`: generated calculation results.
- `图1_匀速巡航_坡度.png`, `图2_匀速巡航_逆风.png`, and `图3_加速与巡航.png`: publication figures.
- `阅读说明.txt`: bundle contents and interpretation limits.

## Build, Test, and Development Commands

Use Python 3.10 or newer; only the standard library is required. No package installation or build step is configured.

```powershell
python .\复算模型.py
```

This regenerates both CSV files and the parameter JSON beside the script, runs numerical assertions, and prints scenario summaries. Change defaults in the `Model` dataclass; the script does not read settings from JSON. The command does not regenerate PNG figures.

Preview Markdown in an editor supporting LaTeX, tables, and relative images before publication.

## Coding Style & Naming Conventions

Use four-space indentation, type hints, `snake_case` functions and variables, and `PascalCase` classes. Follow the existing standard-library approach; no formatter or linter is configured. Keep units explicit in names, such as `_kmh`, `_kg`, and `_w`. Preserve descriptive Chinese artifact names, UTF-8 text, and UTF-8-with-BOM CSV exports.

## Testing Guidelines

There is no separate test framework or coverage threshold. Run the model without Python's `-O` flag so assertions execute. Existing checks verify power-to-speed inversion and compare analytical acceleration power against numerical integration. Expect `All numerical checks passed.` on success.

For model changes, add relevant assertions and reconcile affected article tables, appendix values, data exports, and figures. Numerical checks establish implementation consistency, not experimental validation.

## Commit & Pull Request Guidelines

No Git history is available in this directory, so no established commit convention can be inferred. Use concise imperative summaries, such as `Correct acceleration power calculation`.

Pull requests should explain the changed assumptions or equations, identify affected artifacts, and report validation results. Include figure previews for visual changes and link related issues when applicable. Preserve reference numbering and distinguish crank mechanical power, total-mass W/kg, and rider-mass W/kg.
