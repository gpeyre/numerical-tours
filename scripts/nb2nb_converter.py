import json
import re
import os
import sys

import nt_conversion_lib as lib


class Converter(object):

    def __init__(self, fname):
        with open(fname) as fid:
            self.tree = json.load(fid)
        self.fname = os.path.basename(fname)
        self.name, ext = os.path.splitext(self.fname)
        if 'matlab' in fname:
            self.dest_type = 'python'
        else:
            self.dest_type = 'matlab'

    def convert(self, destination_dir=None):
        # walk the tree, removing output and replacing the code
        # also replace the first code block and the installation block
        if not destination_dir:
            if self.dest_type == 'matlab':
                destination_dir = '../matlab'
            else:
                destination_dir = '../python'

        ws = self.tree['worksheets'][0]['cells']
        if self.dest_type == 'python':
            intro_func = self.get_python_intro
            trans_func = mat2py
        else:
            intro_func = self.get_matlab_intro
            trans_func = py2mat

        first_code = True
        for item in ws:
            if item['cell_type'] == 'code':
                item['outputs'] = []
                if first_code:
                    item['input'] = intro_func()
                    first_code = False
                else:
                    source = [trans_func(line) for line in item['input']]
                    if self.dest_type == 'python' and source[0] == '%%matlab':
                        source = source[1:]
                    if self.dest_type == 'matlab':
                        source.insert(0, '%%matlab')
                    item['input'] = source

            elif item['cell_type'] == 'markdown':
                if item['source'][0].startswith(
                        "*Important:* Please read the [installation page]"):
                    item['source'][0] = lib.INTRO % self.dest_type

        path = os.path.join(destination_dir, self.fname)
        with open(path, 'w') as fid:
            json.dump(self.tree, fid, indent=2, sort_keys=True)

    def get_python_intro(self):
        setup = [
            'from __future__ import division',
            'import nt_toolbox as nt',
            'from nt_solutions import %s as solutions' % self.name,
            '%matplotlib inline',
            '%load_ext autoreload',
            '%autoreload 2']

        return setup

    def get_matlab_intro(self):
        setup = ['%load_ext pymatbridge  # <- put me in my own cell',
                 '%%matlab',
                 "addpath('solutions/%s')" % self.name]
        return setup

    @staticmethod
    def _reformat(text):
        lines = text.splitlines()
        return [l.strip() for l in lines]


def mat2py(line):
    if line.rstrip().endswith(';'):
        line = line.rstrip()[:-1]

    for (pattern, repl) in PY_REPLS:
        line = re.sub(pattern, repl, line)

    if line.lstrip().startswith('for '):
        rest = line.lstrip().split('for ')[1]
        var, _, rest = rest.partition('=')
        rng, _, comment = rest.partition('%')
        line = 'for %s in %s:' % (var, rng.rstrip())
        if comment:
            line += '  # %s' % comment

    # remove right side spacing
    for op in ['\\(', '\[', '{']:
        line = re.sub('%s\s+' % op, op[-1], line)

    # remove left side spacing
    for op in ['\\)', '\]', ':', '}']:
        line = re.sub('\s+%s' % op, op[-1], line)

    # add a space on the left
    for op in [r'\+', '<', '>', '=']:
        line = re.sub(r'(\S)%s' % op,
                      lambda m: '%s %s' % (m.groups()[0], op[-1]),
                      line)

    # add a space on the right
    for op in ['/', r'\+', '=', ':', ',', ';']:
        line = re.sub(r'%s(\S)' % op,
                      lambda m: '%s %s' % (op[-1], m.groups()[0]),
                      line)

    line = re.sub('< =', '<=', line)
    line = re.sub('> =', '>=', line)
    return line


def py2mat(line):
    for (pattern, repl) in MAT_REPLS:
        line = re.sub(pattern, repl, line)

    if line.lstrip().startswith('##'):
        line = '%%' + line.lstrip()[2:]
    elif line.lstrip().startswith('#'):
        line = '%' + line.lstrip()[1:]

    return line


if __name__ == '__main__':
    usage = 'nb2nbconverter.py fname [destination_dir]'

    if len(sys.argv) >= 3:
        c = Converter(sys.argv[1])
        c.convert(sys.argv[2])
    elif len(sys.argv) >= 2:
        c = Converter(sys.argv[1])
        c.convert()
    else:
        print(usage)
