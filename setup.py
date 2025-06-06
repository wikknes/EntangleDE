from setuptools import setup, find_packages

setup(
    name="EntangleDE",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "numpy>=1.20.0",
        "pandas>=1.3.0",
        "qiskit>=0.34.0",
        "scikit-learn>=1.0.0", 
        "matplotlib>=3.4.0",
        "seaborn>=0.11.0",
        "scipy>=1.7.0",
        "statsmodels>=0.13.0",
    ],
    entry_points={
        'console_scripts': [
            'qdta=main:main',  # Keep QDTA as the short CLI command
        ],
    },
    author="Your Name",
    author_email="your.email@example.com",
    description="Quantum Differential Gene Expression Analysis along Pseudotime Trajectories",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/wikknes/EntangleDE",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Artificial Intelligence",
    ],
    python_requires='>=3.7',
)
