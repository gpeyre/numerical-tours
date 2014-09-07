import os
import numpy as np
import nb_converter

def convert_all(matlabdir="../../matlab/", categ='audio'):
	"""
		Process all matlab .m file and convert them into notebooks.
	"""
	A = os.listdir(matlabdir)
	for i in np.arange(1,len(A)):
		a = A[i]
		if a[-2:]=='.m' and a.startswith(categ):
			nb_converter.convert(matlabdir + a)
			