import os
import re
import sys

from template import Notebook

SECTION_HEADING = re.compile('\A%% \w')

MATH_REPLS = [(re.compile(r'\\\['), '$$'),  # replace latex delimiters
              (re.compile(r'\\\]'), '$$'),
              (re.compile(r'\\\('), '$'),
              (re.compile(r'\\\)'), '$'),
              ]

CODE_REPLS = [(re.compile('@\((.*?)\)'),  # replace anon funcs with lambdas
               lambda m: 'lambda %s: ' % m.groups()[0]),
              (re.compile('\A\W*?end\W*?\Z'), ''),  # kill "end" lines
              (re.compile('\A\W*?clf\W*?\Z'), ''),  # kill "clf" lines
              ]


def convert(fname):
    """Convert an m file into an IPython notebook

    Takes most of the grunt work out of translating a numerical-tours
    matlab file into an ipython notebook.  It will auto-populate the
    notebook with headings, markdown, code, and excercises as apppropriate.
    Ideally, all one then has to do is translate the code itself.

    Some code transformations are in place (see CODE_REPLS and _parse_code 
    in this file).
    """
    with open(fname) as fid:
        lines = fid.readlines()

    nb = Notebook()
    nb.add_heading(lines[0][3:].rstrip())

    state = 'markdown'
    out_lines = []
    for line in lines[1:]:
        new_state, new_line = parse_line(line, state)
        if not new_state == state:
            get_section(nb, state, out_lines)
            out_lines = [new_line]
            state = new_state
        else:
            out_lines.append(new_line)

    fname = fname.replace('.m', '.ipynb')
    fname = os.path.basename(fname)
    dname = os.path.dirname('__file__')
    path = os.path.join(os.path.abspath(dname), fname)

    nb.save(path)


def parse_line(line, state):
    new_line = line.decode('utf-8', 'ignore')
    new_state = state
    if state == 'excercise':
        if new_line.startswith('%'):
            new_state = 'excercise'
            new_line = _parse_markdown(new_line)
        elif new_line.startswith('%EXO'):
            new_state = 'markdown'
            new_line = ''
        else:
            new_state = 'excercise'
            new_line = ''
    elif re.match(SECTION_HEADING, new_line):
        new_state = 'section'
        new_line = new_line[3:]
    elif new_line.startswith('%EXO'):
        new_state = 'excercise'
        new_line = ''
    elif new_line.startswith('%'):
        new_state = 'markdown'
    elif not new_line.strip():
        new_state = 'markdown'
        new_line = ''
    else:
        new_state = 'code'
    if new_state == 'markdown':
        new_line = _parse_markdown(new_line)
    elif new_state == 'code':
        new_line = _parse_code(new_line)
    return new_state, new_line


def _parse_markdown(line):
    while line.startswith('%'):
        line = line[1:]
    if line.startswith(' '):
        line = line[1:]
    if '\\' in line:
        for (pattern, repl) in MATH_REPLS:
            line = re.sub(pattern, repl, line)
    return line.rstrip()


def _parse_code(line):
    if line.rstrip().endswith(';'):
        line = line.rstrip()[:-1]
    for (pattern, repl) in CODE_REPLS:
        line = re.sub(pattern, repl, line)
    if line.lstrip().startswith('for '):
        rest = line.lstrip().split('for ')[1]
        var, _, rest = rest.partition('=')
        rng, _, comment = rest.partition('%')
        line = 'for %s in %s:' % (var, rng.rstrip())
        if comment:
            line += '  # %s' % comment
    return line


def get_section(nb, state, out_lines):
    if state == 'section':
        nb.add_heading(out_lines, level=2)
    elif state == 'excercise':
        out_lines = [_parse_markdown(l) for l in out_lines]
        nb.add_exercise(out_lines)
    elif state == 'markdown':
        if len(out_lines) == 1 and not out_lines[0]:
            return
        nb.add_markdown(out_lines)
    elif state == 'code':
        nb.add_code(out_lines)

if __name__ == '__main__':
    if len(sys.argv) >= 2:
        fname = sys.argv[1]
    else:
        fname = '../matlab/meshwav_1_subdivision_curves.m'
    convert(fname)
