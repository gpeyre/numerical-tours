import json
import re


def latex_fixer(obj):
    obj = obj.groups()
    # matches single dollar sign
    if obj[3]:
        return r'\\\(' + obj[4] + r'\\\)'
    # matches double dollar sign
    else:
        return r'\\\[' + obj[1] + r'\\\]'


MATH_REPL = re.compile('(\$\$)(.*?)(\$\$)|(\$)(.*?)(\$)')


def convert(nbfile, mfile=None):
    """Convert a Numerical Tours Matlab notebook file to an m-file.
    """

    output = ''
    with open(nbfile, 'rb') as fid:
        data = json.load(fid)

    ws = data['worksheets'][0]
    for cell in ws['cells']:
        source = cell.get('source', cell.get('input', ''))
        if cell['cell_type'] == 'heading':
            if source[0].startswith('Excercise :'):
                output += '%EX0\n'
            else:
                output += '%% %s' % source[0]
            output += '%'.join(source[1:])
        elif cell['cell_type'] == 'markdown':
            if source[0] == '## Insert your code here.':
                source = ['%EX0']
            source = [line for line in source
                        if not r'\newcommand' in line]
            source = [re.sub(MATH_REPL, latex_fixer, line) for line in source]
            output += '% ' + '% '.join(source)
        elif cell['cell_type'] == 'code':
            output += ''.join(source)
        output += '\n\n'

    if mfile is None:
        mfile = nbfile.replace('.ipynb', '.m')

    with open(mfile, 'wb') as fid:
        fid.write(output)


if __name__ == '__main__':
    convert("../python/denoisingsimp_2b_linear_image.ipynb")

