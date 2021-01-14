from distutils.core import setup

if __name__== '__main__':
    setup(include_package_data=True,
          description='Python OMPA',
          long_description="""Python package for conducting OMPA (Optimum Multiparameter Water Mass Analysis)""",
          url='https://github.com/nitrogenlab/pyompa',
          version='0.1.0.0',
          packages=['pyompa'],
          setup_requires=[],
          install_requires=['numpy', 'pandas', 'cvxpy', 'scipy'],
          extras_require={'altair': ['altair']},
          scripts=[],
          name='pyompa')
