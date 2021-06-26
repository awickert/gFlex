#! /usr/bin/env python3

# modified from:
# https://coderwall.com/p/qawuyq/use-markdown-readme-s-in-python-modules

import pypandoc
import os

output = pypandoc.convert_file('README.md', 'rst')

f = open('README.txt','w+')
f.write(output)
f.close()
