import os
import m2nb_converter

def convert_all(matlabdir="../matlab/m_files/", ntype='python',
                out_dir='../python'):
    """
    Process all matlab .m file and convert them into notebooks.
    """
    files = os.listdir(matlabdir)
    for fname in files:
        if fname.endswith('.m'):
            path = os.path.join(matlabdir, fname)
            converter = m2nb_converter.Converter(path)
            converter.convert(out_dir, ntype)


if __name__ == '__main__':
    convert_all(out_dir='../python/todo')
    #convert_all(ntype='matlab', out_dir='../matlab/todo')
