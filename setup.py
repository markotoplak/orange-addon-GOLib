from distutils.core import setup, Extension

module=Extension("_GOLib", sources=["go.c"])

setup(name="GOLib",
      version="0.1",
      description="A library for handling gene ontologies",
      author="Ales Erjavec",
      author_email="ales.erjavec@fri.uni-lj.si",
      ext_modules=[module], py_modules=["go"]
      )