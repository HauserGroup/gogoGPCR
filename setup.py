from setuptools import setup, find_packages
import pathlib

here = pathlib.Path(__file__).parent.resolve()

long_description = (here / "README.md").read_text(encoding="utf-8")

setup(
    name="gogoGPCR",
    version="0.1.0",
    description="Association testing in the UK Biobank with Hail and regenie",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/jsture/gogoGPCR",
    author="Jakob S. Madsen",
    author_email="jsmadsen@sund.ku.dk",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.6.5",
        "Programming Language :: Python :: 3 :: Only",
    ],
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    python_requires=">=3.6,",
    install_requires=[
        "hail==0.2.78",
        "tomli",
        "jupyterlab",
    ],
    extras_require={
        "dev": [
            "black",
            "isort",
        ],
    },
    # entry_points={  # Optional
    #     "console_scripts": [
    #         "sample=sample:main",
    #     ],
    # },
    project_urls={
        "Bug Reports": "https://github.com/jsture/gogoGPCR/issues",
        "DNAnexus environment": (
            "https://ukbiobank.dnanexus.com/panx/tools/workers/jupyterLab"
        ),
        "Source": "https://github.com/jsture/gogoGPCR",
    },
)