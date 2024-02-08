from setuptools import find_packages
from distutils.core import setup

if __name__== '__main__':
    setup(include_package_data=True,
          description='Python OMPA',
          long_description="""Python package for conducting OMPA (Optimum Multiparameter Water Mass Analysis)""",
          author="Avanti Shrikumar",
          author_email="avanti.shrikumar@gmail.com",
          url='https://github.com/nitrogenlab/pyompa',
          version='0.4.2.1',
          packages=find_packages(),
          setup_requires=[],
          install_requires=['numpy', 'pandas', 'cvxpy', 'scipy', 'toml'],
          extras_require={'altair': ['altair']},
          scripts=['scripts/run_ompa_given_config'],
          name='pyompa')
