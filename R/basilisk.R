#' Python environments
autogenes <- BasiliskEnvironment(
    envname = "autogenes",
    pkgname = "StatescopeR",
    packages = c(
        "autogenes==1.0.4",
        "python==3.7", # https://github.com/theislab/AutoGeneS/issues/27
        "pandas==1.3.5", "anndata==0.8.0",
        "numpy==1.21.6", "dill==0.3.7",
        "deap==1.4.3", "scipy==1.7.3",
        "cachetools==5.5.2",
        "scikit-learn==1.0.2",
        "matplotlib==3.5.3"
    ),
    channels = c("anaconda", "bioconda", "conda-forge")
)

deconvolution <- BasiliskEnvironment(
    envname = "deconvolution",
    pkgname = "StatescopeR",
    packages = c(
        "python==3.11", "numba==0.62.1",
        "pandas==1.5.3", "joblib==1.5.2",
        "numpy==1.26.4",
        "scipy==1.16.3",
        "scikit-learn==1.5.2",
        "matplotlib==3.6.3",
        "mkl==2025.3.0",
        "torch==2.9.1", "dill==0.3.4"
    ),
    channels = c("bioconda", "conda-forge", "pytorch")
)

statescope <- BasiliskEnvironment(
    envname = "statescope",
    pkgname = "StatescopeR",
    packages = c(
        "python==3.11",
        "numba==0.62.1",
        "pandas==1.5.3",
        "joblib==1.5.2",
        "numpy==1.23.5",
        "scipy==1.15.3",
        "matplotlib==3.6.3",
        "seaborn==0.13.2",
        "scikit-learn==1.5.2"
    ),
    channels = c("bioconda", "conda-forge")
)
