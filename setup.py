from setuptools import setup

setup(
    name="DBSpace",
    version="0.1.1",
    author="Vineet Tiruvadi",
    author_email="virati@gmail.com",
    packages=["dbspace"],
    package_dir={"": "src"},
    url="http://pypi.python.org/pypi/DBSpace/",
    license="LICENSE.txt",
    description="Multimodal, Model-based DBS Analysis",
    long_description=open("README.md").read(),
    install_requires=[
        "python==3.8.12",
        "numpy",
    ],
)
