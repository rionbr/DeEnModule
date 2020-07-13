from setuptools import setup, find_packages
from deenmodule import __package__, __description__, __version__, __author__


def readme():
    with open('README.md') as f:
        return f.read()


setup(
    name=__package__,
    version=__version__,
    description=__description__,
    long_description=__description__,
    classifiers=[
        'Development Status :: 4 - Beta',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2.7',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Information Analysis',
    ],
    keywords="entropy pca principal component analysis networks computational biology social science",
    url="http://github.com/rionbr/DeEnModule",
    author=__author__,
    author_email="rionbr@gmail.com",
    license="MIT",
    packages=find_packages(),
    package_data={
        'datasets': [
            'deenmodule.datasets/*.txt',
        ],
    },
    install_requires=[
        'numpy',
        'scipy',
        'pandas',
        'networkx'
    ],
    include_package_data=True,
    zip_safe=False)
