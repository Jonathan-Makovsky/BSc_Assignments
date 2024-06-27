from setuptools import setup, find_packages, Extension

setup(
    name='mykmneassp',
    version='0.1.0',
    author="Yonatan",
    author_email="Yonatan@example.com",
    description="mykmeanssp",
    install_requires=['invoke'],
    packages=find_packages(),

    license ='GPL-2',

    classifiers =[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: GNU General Public License v2 (GPLv2)',
        'Natural language :: English',
        'Programming Language :: python :: 3 :: Only',
        'Programming Language :: Python :: Implementation :: CPython',
    ],
    ext_modules=[
        Extension(
            'mykmeanssp',
            ['spkmeansmodule.c', 'spkmeans.c'],
        ),
    ]
)