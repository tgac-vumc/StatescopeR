library(basilisk)
autogenes <- BasiliskEnvironment(envname="autogenes",
                              pkgname="StatescopeR",
                              packages=c('python==3.10.4', "autogenes==1.0.4",
                                         'pandas==1.5.3','anndata==0.9.2',
                                         'numpy==1.24.4', 'dill==0.3.8',
                                         'deap==1.4.1', 'scipy==1.14.1',
                                         'cachetools==4.2.2', 'scikit-learn==1.1.3',
                                         'matplotlib==3.9.3'),
                              channels = c('bioconda', 'conda-forge')
)
