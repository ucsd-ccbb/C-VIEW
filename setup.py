from setuptools import setup, find_packages

setup(
    name="sarscov2-sequencing-pipeline",
    author="UCSD Center for Computational Biology and Bioinformatics",
    author_email="ccbb@health.ucsd.edu ",
    packages=find_packages(),
    version='0.4.0',
    url="https://github.com/ucsd-ccbb/covid_sequencing_analysis_pipeline",
    description="A high-throughput data processing pipeline for SARS-CoV-2 variant identification",
    license='MIT License',
    package_data={'covid_sequencing_analysis_pipeline':
                  [
                      'pipeline/*.*',
                      'qc/*.*',
                      'reference_files/*.*',
                      'tests/data/*.*'
                  ]}
)
