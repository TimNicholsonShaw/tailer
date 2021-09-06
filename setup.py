from setuptools import setup, find_packages
import pathlib

here = pathlib.Path(__file__).parent.resolve()

long_description = (here / 'README.md').read_text(encoding='utf-8')

setup(
    name='jla-tailer', 
    version='0.0.11',
    description='Tool to find 3\' tailing of non-coding RNAs',
    long_description=long_description,
    long_description_content_type='text/markdown', 
    url='https://github.com/TimNicholsonShaw/tailer',  
    author='Tim Nicholson-Shaw',
    author_email='timnicholsonshaw@gmail.com',
    classifiers=[ 
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3 :: Only',
    ],
    keywords='Sequencing, ncRNAs, Bioinformatics',
    package_dir={'': 'src'},
    packages=find_packages(where='src'),
    python_requires='>=3.6',
    install_requires=['pysam', 'Bio', 'gffutils'],

    entry_points={
        'console_scripts': [
            'Tailer=Tailer.Tailer:main',
        ],
    },

    project_urls={ 
        'Lab Website': 'https://labs.biology.ucsd.edu/lykkeandersen/index.html',

    },
)