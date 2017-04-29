from setuptools import setup
from setuptools import find_packages

setup(name='mesh2mesh',
      version='0.9.0',
      description='Mesh-to-Mesh Field Projection for Finite Element Analysis',
      url='mesh2mesh.com',
      download_url='https://github.com/zifu-wang/mesh2mesh',
      author='Zifu Wang',
      author_email='z@mesh2mesh.com',
      license='GNU General Public License (GPL)',
      platforms = ["Linux", "Mac OS-X", "Unix"],
      packages=find_packages(exclude=['build','test']),
      package_data={'': ['*.tab', '*.so']},
      install_requires=[],
      )

