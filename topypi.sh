#! /usr/bin/env sh

python3 setup.py sdist bdist_wheel
twine upload dist/*
