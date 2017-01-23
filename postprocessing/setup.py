from distutils.core import setup

setup(name='pypp',
      version='0.0',
      description='post-processing eigendecompositions computed with paladin',
      url='https://github.com/michael-a-hansen/paladin',
      author='Mike Hansen',
      author_email='mike.hansen.utah@gmail.com',
      license='MIT',
      packages=['pypp'],
      install_requires=['numpy',
                        'scipy',
                        'matplotlib'],
      zip_safe=False)
