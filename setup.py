from setuptools import setup, find_packages

with open('README.md') as f:
    long_description = f.read()

setup(name="protclus",
    version="0.7.0",
    description="protclus: Library for protein complex discovery on PPI networks",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/paulmorio/protclus",
    author="Paul Scherer",
    author_email="paul.geo2dr@gmail.com",
    license="MIT",
    install_requires=['numpy', 'tqdm', 'networkx', 'py27hash'],
    
    test_suite="nose.collector",
    tests_require=["pytest", "pytest-cov", "mock"],
    packages=find_packages(),
    include_package_data=True,
    classifiers=['Topic :: Scientific/Engineering :: Artificial Intelligence',
                 'License :: OSI Approved :: MIT License',
                 'Natural Language :: English',
                 'Programming Language :: Python :: 3.6',
                 'Operating System :: OS Independent'],
    zip_safe=False
    )