from setuptools import setup, find_packages

setup(
    name='pycircdb',
    version='0.1.0',
    install_requires=[
        'click',
        'rich_click',
        'pyarrow',
        'pandas',
        'dask',
        'dask[distributed]',
        'networkx'
    ],
    entry_points={
        "console_scripts": [
            "pycircdb=pycircdb.__main__:run_pycircdb",
        ],
    },
)
