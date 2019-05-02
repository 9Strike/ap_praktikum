from distutils.core import setup, Extension

sigval = Extension('sigval', sources=['src/sigval.c'], include_dirs=['pyinc'])
setup(
  name = 'Significant value',
  version='1.0',
  description='Parses a value and its error to a string with only the significant part',
  ext_modules = [sigval]
)
