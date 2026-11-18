# Local software verification

Verified on 18 November 2026 on Windows with Python 3.13.16,
NumPy 2.5.4 and base R 4.6.1.

`python -B run_tests.py`: **11 unit tests passed**, zero failures/errors.
All 2,000 Python/base-R replicate differences and both percentile endpoints agreed within 1e-12; the independent exact-discordance hand example was checked in both implementations.
Python source parsed successfully. Staged source hashes were unchanged by execution.
Only reviewed invented inputs and temporary outputs were used; no clinical input,
source-data profiling, model training or dependency installation was performed.

These are same-machine software/arithmetic checks, not clinical validation or an
independent human reproduction. An independent language implementation is not an
independent clinical reviewer. The preserved clinical HOLD/FALSE gates remain unchanged.
No hosted-CI, external-platform test or clean installation is claimed.

The positive RELEASE_MANIFEST.json identifies the exact intended release payload.
Source provenance hashes distinguish existing modules from new documentation/runners.
