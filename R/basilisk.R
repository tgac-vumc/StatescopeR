#' Python environments
statescope <- BasiliskEnvironment(
    envname = "statescope",
    pkgname = "StatescopeR",
    packages = c(
        "python==3.11",
        "anndata==0.11.4",
        "deap==1.4.3",
        "cachetools==7.0.5",
        "numba==0.62.1",
        "pandas==1.5.3",
        "joblib==1.5.3",
        "numpy==1.23.5",
        "scipy==1.15.3",
        "matplotlib==3.6.3",
        "seaborn==0.13.2",
        "scikit-learn==1.5.2",
        "torch==2.8.0", "dill==0.3.4",
        "statescope-autogenes==1.0.4.post2",
        "psutil==7.2.2"
    ),
    channels = c("anaconda", "bioconda", "conda-forge", "pytorch")
)
