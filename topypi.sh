#! /usr/bin/env sh

python register.py
python setup.py register -r pypi
rm README.txt
python setup.py sdist upload -r pypi
