## StatescopeR 0.99.0 (3-6-2025)

-   Initial Bioconductor submission

## StatescopeR 0.99.1-0.99.6 (16-7-2025 - 26-7-2025)

-   Fixing Basilisk issues

## StatescopeR 0.99.7-0.99.21 (28-7-2025 - 23-7-2025)

-   Adjusting Vignette & Examples to make Check process quick enough to not
timeout

## StatescopeR 0.99.22 (22-10-2025)

- Revision after review, in short: Some documentation changes, removal of custom
functions in favor of existing good implementations & removal of new classes in 
favor of adding to SummarizedExperiment metadata, see 
https://github.com/Bioconductor/Contributions/issues/3838#issuecomment-3261548846
for more details.

## StatescopeR 0.99.23 & 0.99.24 (23-10-2025)
-   Minimizing Vignette & Examples to make Check process quick enough to not
timeout


## StatescopeR 0.99.25-0.99.31 (26-11-2025)
-   R dependency >= 4.6.0 reflecting new Bioconductor release and minimized
number of samples for deconvolution for quicker vignette/examples

## StatescopeR 0.99.32 (10-12-2025)
-   Adjusted to use BLADE code with pytorch (quicker version: https://github.com/tgac-vumc/Statescope/blob/master/src/BLADE_Deconvolution/BLADE.py)
and made vignette less minimal

## StatescopeR 0.99.33 (24-12-2025)
-   Added function to fetch signatures from https://github.com/tgac-vumc/StatescopeData and did some code cleaning

## StatescopeR 0.99.34 (28-1-2026)
-   Adjusted vignette based on supervisor comments

## StatescopeR 0.99.35 (1-4-2026)
-   Limited basilisk environments using Statescope-autogenes

## StatescopeR 0.99.36 (8-4-2026)
-   Removed mkl dependency for mac support and removed .bbsoptions for upcoming R universe build system
