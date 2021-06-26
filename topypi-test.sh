#! /usr/bin/env sh

twine upload --repository-url https://test.pypi.org/legacy/ dist/*
#python setup.py sdist upload -r pypitest
