# Verification

Local release candidate v0.4.0 was checked on 25 September 2026.

## Candidate checks completed in this environment

- Package metadata reports version `0.4.0`.
- The complete bundled test suite passed with `PYTHONPATH=src python -m pytest -q`: **13 tests passed**.
- The v0.4.0 revision bundle includes the collected BJP revision outputs, exact run specification, final Hauser et al. functional-repertoire annotation, G12/13 census, 92-feature dictionary, 408-comparison table, four-fold coverage sensitivity files, and revised manuscript figure assets.
- SHA-256 checksums and a full file manifest are included.

## Environment limitation

A normal editable `pip install -e '.[test]'` could not be completed in this offline container because pip attempted to download the declared build dependency `setuptools>=69`. This is an environment/network limitation, not a test failure. Re-run the standard installation and GitHub Actions matrix after pushing the candidate to GitHub.

## Scientific/reproduction boundary

The specimen-level audit distinguishes coordinate-level receptor-only models from experimentally partner-free specimens. The compact revised SI Tables S22-S26 are included. The full raw specimen-provenance record (complete deposited titles, raw EM sample descriptions, particle/shared-particle evidence and exact specimen-audit script) should be added to the Zenodo/source-data release if it is not already present in the author's working archive.
