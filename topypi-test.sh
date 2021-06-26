#! /usr/bin/env sh

python3 register.py
rm README.txt
twine upload --repository-url https://test.pypi.org/legacy/ dist/*
#python setup.py sdist upload -r pypitest
