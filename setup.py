from distutils.core import setup

setup(
    name='Oxygen_Project',
    version='0.1.0',
    author='Abigail Caron',
    author_email='abigailc@mit.edu',
    packages=['oxy_mod'],
    scripts=['bin/Overall_DTL_Detector.py'],
    url='http://pypi.python.org/pypi/Oxygen_Project/',
    license='LICENSE.txt',
    description='Useful towel-related stuff.',
    long_description=open('README.txt').read(),
    install_requires=[
        "os",
        "caldav == 0.1.4",
    ],
)