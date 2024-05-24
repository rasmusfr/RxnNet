from setuptools import setup, find_packages


with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='rxnnet',
    long_description=readme,
    author='Rasmus Fromsejer',
    url='https://github.com/rasmusfr/RxnNet',
    license=license,
    packages=find_packages(),
    install_requires=['pandas~=2.2.0', 'numpy~=1.26.4', 'wquantiles~=0.6', 'chempy~=0.8.3', 'pymatgen~=2024.3.1'],
    version='0.0.1'
)