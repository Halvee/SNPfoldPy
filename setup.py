from distutils.core import setup
from os.path import join as pjoin



setup(name='SNPfoldPy',
      version='1.02',
      description='tool for quantifying the predicted level of global secondary structure change brought about by a SNP or SNPs in an RNA sequence.',
      author='Matt Halvorsen',
      author_email='mhalvors1@gmail.com',
      url='',
      packages=['SNPfold','SNPfold.BasicFunctions'],
      scripts=[pjoin('bin', 'SNPfold_commandline.py')]
      )

