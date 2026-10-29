# module-3-assignment

MPCS 56430 Module 3 assignment. The instructions are in the notebooks in
`notebooks/`.

## Setup

The project uses [uv](https://docs.astral.sh/uv/), as described in the course
handbook's *Python projects* chapter. From the project folder:

```
uv sync
uv run pytest
```

`uv sync` creates `.venv`, installs the Python version in `.python-version`,
the dependencies at the versions recorded in `uv.lock`, and this project's
package, `module_3_assignment`. The tests should pass before you change anything.

In VS Code, open a notebook, choose **Select Kernel**, then **Python
Environments**, and pick the interpreter in `.venv`.

## Commands

| Task | Command |
|---|---|
| Run all tests | `uv run pytest` |
| Run one test file | `uv run pytest tests/test_package.py -q` |
| Lint | `uv run ruff check src tests` |
| Fix what the linter can fix | `uv run ruff check --fix src tests` |
| Add a package | `uv add <package>` |
| Build API docs into `docs/` | `uv run pdoc src/module_3_assignment -o docs` |

Commit `pyproject.toml` and `uv.lock` whenever you add a package, so the
environment can be rebuilt with `uv sync --locked`.

## Layout

| Path | Contents |
|---|---|
| `src/module_3_assignment/` | Your reusable code, imported as `module_3_assignment` |
| `tests/` | Tests for that code |
| `notebooks/` | Assignment notebooks; they import from `module_3_assignment` |
| `data/` | Input data |
| `docs/` | Generated documentation |
| `CLAUDE.md` | The course's AI policy, for Claude Code |

## Continuous integration

`.github/workflows/ci.yml` runs on every push to `main` and every pull
request. It rebuilds the environment with `uv sync --locked`, then runs
`ruff` and `pytest`. A red check on GitHub means one of those failed.
