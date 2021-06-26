#! /usr/bin/env sh

python setup.py sdist bdist_wheel
twine upload -r testpypi dist/*
