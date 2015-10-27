#! /usr/bin/env sh

python register.py
python setup.py register -r pypitest
rm README.txt
python setup.py sdist upload -r pypitest
