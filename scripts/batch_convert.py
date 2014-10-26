import os
from . import nb_converter


def convert_all(matlabdir="../matlab/", ntype='python', 
               out_dir='../python'):
    """
    Process all matlab .m file and convert them into notebooks.
    """
    files = os.listdir(matlabdir)
    for fname in files:
        if fname.endswith('.m'):
            nb_converter.convert(matlabdir + fname, ntype)


if __name__ == '__main__':
    convert_all()
    convert_all(ntype='matlab', out_dir='../matlab/notebooks')
