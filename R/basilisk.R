#' Python environments
autogenes <- BasiliskEnvironment(envname="autogenes",
                              pkgname="StatescopeR",
                              packages=c('python==3.10.4', "autogenes==1.0.4",
                                         'pandas==1.5.3','anndata==0.9.2',
                                         'numpy==1.24.4', 'dill==0.3.8',
                                         'deap==1.4.1', 'scipy==1.14.1',
                                         'cachetools==4.2.2',
                                         'scikit-learn==1.1.3',
                                         'matplotlib==3.9.3'),
                              channels = c('bioconda', 'conda-forge')
)

deconvolution <- BasiliskEnvironment(envname="deconvolution",
                                 pkgname="StatescopeR",
                                 packages=c('python==3.11', 'numba==0.59.1',
                                            'pandas==1.5.3', 'joblib==1.4.2',
                                            'numpy==1.23.5',
                                            'scipy==1.14.1',
                                            'scikit-learn==1.1.3',
                                            'matplotlib==3.9.3'),
                                 channels = c('bioconda', 'conda-forge')
)

statescope <- BasiliskEnvironment(envname="statescope",
                                     pkgname="StatescopeR",
                                     packages=c('python==3.11',
                                                'numba==0.59.1',
                                                'pandas==1.5.3',
                                                'joblib==1.4.2',
                                                'numpy==1.23.5',
                                                'scipy==1.14.1',
                                                'matplotlib==3.9.3',
                                                'cvxopt==1.3.2',
                                                'seaborn==0.13.2',
                                             'scikit-learn==1.6.0'),
                                     channels = c('bioconda', 'conda-forge')
)

