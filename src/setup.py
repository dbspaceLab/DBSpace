from setuptools import setup

setup(
    name='DBSpace',
    version='0.9.0',
    author='Vineet Tiruvadi',
    author_email='virati@gmail.com',
    packages=['DBSpace', 'DBSpace.test'],
    url='http://pypi.python.org/pypi/DBSpace/',
    license='LICENSE.txt',
    description='Multimodal, Model-based DBS Analysis',
    long_description=open('README.md').read(),
    install_requires=[
        "numpy",
        "pytest",
    ],
)
