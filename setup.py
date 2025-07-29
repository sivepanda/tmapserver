from setuptools import setup, find_packages
print(find_packages())

setup(
    name='tmapserver',
    version='0.0.2',
    packages=find_packages(include=['tmapserver']),
    install_requires=[
        'numpy',
        'tqdm',
        'STHD',
        'pandas',
        'anndata',
        'scipy',
        'subprocess-tee',  # For better subprocess handling
        'psutil',          # For process management
        'requests',        # For checking server status
    ],
    entry_points={
        'console_scripts': [
            'tmapserver=tmapserver.main:main',
        ],
    },
)
