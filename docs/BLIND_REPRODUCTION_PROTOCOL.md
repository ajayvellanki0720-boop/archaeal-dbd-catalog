# Blind reproduction protocol

An analyst who did not develop the pipeline should receive only the archived release, the manuscript, and this protocol.

1. Record operating system, hardware, container runtime, and start time.
2. Verify `MANIFEST.sha256` before running anything.
3. Follow the README without contacting the authors unless progress is impossible.
4. Run the smoke test, automated tests, release audit, and every analysis exposed by `archaeal-dbd all`.
5. Compare every regenerated table byte-for-byte where deterministic and numerically within documented tolerance otherwise.
6. Record every undocumented dependency, ambiguous instruction, manual intervention, warning, runtime, and failed command.
7. Classify discrepancies as documentation, environment, nondeterminism, missing data, or analytical mismatch.
8. Deposit the signed report and machine-readable comparison table with the release.

Success requires all tests and invariants to pass and no unexplained numerical discrepancy in a manuscript-relevant result. This protocol is ready, but independence requires a person outside the analysis team and therefore cannot be self-certified here.
