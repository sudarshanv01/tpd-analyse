import json
from setuptools import setup, find_packages

def setup_package():

    _filename_setup_json = 'setup.json'
    _filename_description = 'README.md'

    with open(_filename_setup_json, 'r') as handle:
        setup_json = json.load(handle)

    with open(_filename_description, 'r') as handle:
        description = handle.read()

    setup(
        **setup_json
    )


if __name__ == '__main__':
    setup_package()