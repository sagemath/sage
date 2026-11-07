# Homebrew dependency bottle CI

The `Homebrew bottles` workflow builds the checked-out Sage sources against
`dimpase/tap/sagemath-deps` and `dimpase/tap/sagemath-docbuild-deps`. It runs on
the self-hosted ARM64 macOS runner labeled `darwin-27`, using Homebrew's Python
3.14 and at most four build or test workers. The dependency formulae provide
Python and native libraries, including NumPy and SciPy built with Accelerate.
The Sage build disables documentation and Meson subproject fallbacks.

Pull requests and pushes to `develop` or `master` run all pytest tests under
`src/sage` and a smoke selection of doctests covering arithmetic, integers,
double matrices, PARI, eclib, GAP, and Maxima. Test failures fail the job.
Import checks verify that both Python and compiled Sage modules come from the
checkout and its build directory, and that dependency modules come from the
Homebrew bundles. Each run creates and removes its own venv and build directory;
the dependency installation and compiler cache persist between runs. Jobs in
this workflow are serialized because they share the Homebrew installation.

The workflow also supports manual dispatch. For example, to test the head of
an upstream PR:

```sh
gh workflow run homebrew-bottles.yml --repo sagemath/sage \
  -f sage_repository=sagemath/sage -f sage_ref=refs/pull/12345/head
```

Leave `sage_repository` and `sage_ref` empty to test the workflow's own commit.
Set `suite=all` to run all doctests under `src/sage`, or supply space-separated
`doctest_paths` to select existing paths under `src/sage`. Explicit paths take
precedence over `suite`; pytest tests always run. As in Sage's other workflows,
open PRs labeled `p: CI Fix` are applied before building.

Every run uploads an artifact with dependency versions, tap commits, data
locations, import provenance, build and test logs, pytest XML, doctest statistics,
and Meson diagnostics. Artifacts are retained for 14 days. The GitHub job summary
records the tested Sage commit and test selection.
