import os
import re
from setuptools import find_packages

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

version = '1.0.0'

_REQUIREMENTS_FILE = 'REQUIREMENTS.txt'
_README = 'README.md'

def _get_local_file(file_name):
    return os.path.join(os.path.dirname(__file__), file_name)

def _get_description(file_name):
    with open(file_name, 'r') as f:
        _long_description = f.read()
    return _long_description

def _get_requirements(file_name):
    with open(file_name, 'r') as f:
        lines = f.readlines()
    rx = re.compile('^[A-z]')
    requirements = [l for l in lines if rx.match(l) is not None]
    return requirements

setup(
    name='FALCON-polish',
    version=version,
    package_dir={'': '.'},
    packages=[
        'falcon_polish',
        #'falcon_polish.cli',
        'falcon_polish.mains',
        'falcon_polish.pypeflow', # might be moved to FALCON someday
    ],
    license='See LICENSE file.',
    author='PacBio',
    author_email='cdunn@pacificbiosciences.com',
    description='FALCON using pbsmrtpipe.',
    #setup_requires=['nose>=1.0'],
    # Maybe the pbtools-* should really be done in a subparser style
    entry_points={'console_scripts': [
        #'pb-falcon = pbfalcon.run:main',
        #'pb-ini2xml = pbfalcon.ini2xml:main',
    ]},
    install_requires=_get_requirements(_get_local_file(_REQUIREMENTS_FILE)),
    tests_require=['nose'],
    long_description=_get_description(_get_local_file(_README)),
    #classifiers=['Development Status :: 4 - Beta'],
    include_package_data=True,
    zip_safe=False,
    #package_data={'pbfalcon': ['reg-tcs/*.json', ]},
)
