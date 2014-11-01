import os
import json
import nt_conversion_lib as lib


def fix_intro(path):
    with open(path) as fid:
        tree = json.load(fid)

    if 'matlab' in path:
        out_type = 'matlab'
    else:
        out_type = 'python'

    ws = tree['worksheets'][0]['cells']
    intro = ws[0]['source']
    ws[0]['source'] = intro[:2] + [lib.INTRO % out_type] + intro[2:]

    tree['worksheets'][0]['cells'] = ws[:4] + ws[5:]

    with open(path, 'w') as fid:
        json.dump(tree, fid)


def convert_all(dirname):
    """
    Process all matlab .m file and convert them into notebooks.
    """
    files = os.listdir(dirname)
    for fname in files:
        if fname.endswith('.ipynb'):
            path = os.path.join(dirname, fname)
            fix_intro(path)


if __name__ == '__main__':
    convert_all('../matlab/')
    #fix_intro('../matlab/sparsity_9_sparsespikes_cbp.ipynb')
