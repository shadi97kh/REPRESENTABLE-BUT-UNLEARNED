# Final audit correction

The first delivery audit incorrectly required float64 subtraction of exported endpoints to equal the original float32 B2 difference within 1e-10. The maximum discrepancy was 1.49011613e-8. All 4,680 neural pair predictions match the original binary32 subtraction exactly after binary32 conversion; CSV parsing contributes less than 1e-16. The corrected check requires exact binary32 equality and a half-ULP real-arithmetic rounding bound plus 1e-15 CSV tolerance. The failed log, original audit source and diagnostic rows remain in the current run. No fit, prediction, metric or scientific conclusion changed.

Delivery subprocess logs now receive a unique attempt directory so retries preserve prior logs. All attempts count as administration, not new fits.
