#!/bin/env python
import os
import sys
# import git
import subprocess

# git.Git("HaMStR").clone("https://github.com/trvinh/HaMStR")
# cmd = "git clone https://github.com/trvinh/HaMStR"
# subprocess.call(cmd, shell = True)
# os.chdir("HaMStR")
subprocess.call("bin/setup.sh")
