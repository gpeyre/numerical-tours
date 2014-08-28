import os
import re
import sys

from nb_template import Notebook

SECTION_HEADING = re.compile('\A%% \w')

MATH_REPLS = [(re.compile(r'\\\['), '$$'),  # replace latex delimiters
              (re.compile(r'\\\]'), '$$'),
              (re.compile(r'\\\('), '$'),
              (re.compile(r'\\\)'), '$'),
              ]

PY_REPLS = [(re.compile('@\((.*?)\)'),  # replace anon funcs with lambdas
             lambda m: 'lambda %s: ' % m.groups()[0]),
            (re.compile('\A\W*?end\W*?\Z'), ''),  # kill "end" lines
            (re.compile('\A\W*?clf\W*?\Z'), ''),  # kill "clf" lines
            ]


class Converter(object):

    """Convert an m file into an IPython notebook

    Takes most of the grunt work out of translating a numerical-tours
    matlab file into an ipython notebook.  It will auto-populate the
    notebook with headings, markdown, code, and exercises as appropriate.
    Ideally, all one then has to do is translate the code itself.

    For python, some code transformations are in place
     (see CODE_REPLS and `parse_code`).
    """

    def __init__(self, fname, ntype='python'):
        self.ntype = ntype.lower()
        name = os.path.basename(fname)
        name = name.replace('.m', '')
        self.nb = Notebook(name, self.ntype)
        self.fname = fname

    def convertt(self):
        with open(self.fname) as fid:
            lines = fid.readlines()

        self.nb.add_heading(lines[0][3:].rstrip())

        state = 'markdown'
        out_lines = []
        for line in lines[1:]:
            new_state, new_line = self.parse_line(line, state)
            if not new_state == state:
                self.get_section(state, out_lines)
                out_lines = [new_line]
                state = new_state
            else:
                out_lines.append(new_line)
        # handle the last section
        self.get_section(state, out_lines)

        fname = self.fname.replace('.m', '.ipynb')
        if self.ntype == 'python':
            fname = os.path.basename(fname)
            dname = os.path.dirname('__file__')
            path = os.path.join(os.path.abspath(dname), fname)
        else:
            dname = os.path.dirname(fname)
            fname = os.path.basename(fname)
            fname = os.path.join(dname, 'notebooks', fname)

        self.nb.save(path)

    def parse_line(self, line, state):
        new_line = line.decode('utf-8', 'ignore')
        new_state = state
        if state == 'excercise':
            if new_line.startswith('%EXO'):
                new_state = 'markdown'
                new_line = ''
            elif new_line.startswith('%'):
                new_state = 'excercise'
                new_line = self.parse_markdown(new_line)
            else:
                new_state = 'excercise'
                new_line = ''
        elif state == 'comment':
            if new_line.startswith('%CMT'):
                new_state = 'markdown'
                new_line = ''
        elif re.match(SECTION_HEADING, new_line):
            new_state = 'section'
            new_line = new_line[3:]
        elif new_line.startswith('%EXO'):
            new_state = 'excercise'
            new_line = ''
        elif new_line.startswith('%CMT'):
            new_state = 'comment'
            new_line = ''
        elif new_line.startswith('%'):
            new_state = 'markdown'
        elif not new_line.strip():
            new_state = 'markdown'
            new_line = ''
        else:
            new_state = 'code'
        if new_state == 'markdown':
            new_line = self.parse_markdown(new_line)
        elif new_state == 'code' and self.ntype == 'python':
            new_line = self.parse_code(new_line)
        return new_state, new_line

    def parse_markdown(self, line):
        while line.startswith('%'):
            line = line[1:]
        if line.startswith(' '):
            line = line[1:]
        if '\\' in line:
            for (pattern, repl) in MATH_REPLS:
                line = re.sub(pattern, repl, line)
        return line.rstrip()

    def parse_code(self, line):
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

    def get_section(self, state, out_lines):
        nb = self.nb
        if state == 'section':
            nb.add_heading(out_lines, level=2)
        elif state == 'excercise':
            out_lines = [self.parse_markdown(l) for l in out_lines]
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
    Converter(fname).convert()
