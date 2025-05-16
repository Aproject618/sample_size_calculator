from setuptools import setup, find_packages

setup(
    name='sample-size-calculator',
    version='0.1.0',
    description='A Python package for calculating required sample sizes for A/B testing using z-tests',
    long_description=open('README.md', encoding='utf-8').read(),
    long_description_content_type='text/markdown',
    author='Ramis Sungatullin',
    author_email='aproject618@gmail.com',
    url='https://github.com/yourusername/sample-size-calculator',
    packages=find_packages(),
    python_requires='>=3.7',
    install_requires=[
        'numpy>=1.18.0',
        'scipy>=1.4.0',
        'plotly>=5.0.0',
    ],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Topic :: Scientific/Engineering :: Information Analysis',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ],
    include_package_data=True,
    license='MIT',
    keywords='ab-testing statistics sample-size z-test',
)
