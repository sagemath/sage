# Homebrew dependency bottle CI

The workflow builds Sage from the checkout with uv and the Homebrew dependency
bottles on macOS and Linux. It runs pytest and the full Sage doctest suite.

To test another Sage ref, dispatch the workflow manually:

```sh
gh workflow run release.yml --repo sagemath/sage -f target=homebrew \
  -f sage_repository=sagemath/sage -f sage_ref=refs/pull/12345/head
```

For a shorter doctest run, pass `-f suite=smoke` or use `-f doctest_paths="..."`
with paths under `src/sage`. The workflow uploads build and test logs as artifacts.
