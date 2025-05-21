#' Python environments
autogenes <- BasiliskEnvironment(
    envname = "autogenes",
    pkgname = "StatescopeR",
    packages = c(
        "autogenes==1.0.4",
        "python==3.6", #https://github.com/theislab/AutoGeneS/issues/27
        "pandas==1.1.5", "anndata==0.7.8",
        "numpy==1.19.5", "dill==0.3.4",
        "deap==1.4.1", "scipy==1.5.4",
        "cachetools==4.2.2",
        "scikit-learn==0.24.2",
        "matplotlib==3.3.4"
    ),
    channels = c('anaconda', "bioconda", "conda-forge")
)

deconvolution <- BasiliskEnvironment(
    envname = "deconvolution",
    pkgname = "StatescopeR",
    packages = c(
        "python==3.11", "numba==0.59.1",
        "pandas==1.5.3", "joblib==1.4.2",
        "numpy==1.23.5",
        "scipy==1.14.1",
        "scikit-learn==1.1.3",
        "matplotlib==3.9.3"
    ),
    channels = c("bioconda", "conda-forge")
)

statescope <- BasiliskEnvironment(
    envname = "statescope",
    pkgname = "StatescopeR",
    packages = c(
        "python==3.11",
        "numba==0.59.1",
        "pandas==1.5.3",
        "joblib==1.4.2",
        "numpy==1.23.5",
        "scipy==1.14.1",
        "matplotlib==3.9.3",
        "cvxopt==1.3.2",
        "seaborn==0.13.2",
        "scikit-learn==1.6.0"
    ),
    channels = c("bioconda", "conda-forge")
)
