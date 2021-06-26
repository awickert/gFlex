#! /usr/bin/env sh

rm -r dist
python3 setup.py sdist bdist_wheel
twine upload dist/*
