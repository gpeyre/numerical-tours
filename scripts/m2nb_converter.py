import os
import re
import sys

from nb_template import Notebook
import nt_conversion_lib as lib
from convert_to_matlab_kernel import convert_to_matlab_kernel


class Converter(object):

    """Convert an m file into an IPython notebook

    Takes most of the grunt work out of translating a numerical-tours
    matlab file into an ipython notebook.  It will auto-populate the
    notebook with headings, markdown, code, and exercises as appropriate.
    Ideally, all one then has to do is translate the code itself.

    For python, some code transformations are in place
     (see CODE_REPLS and `parse_code`).
    """

    def __init__(self, fname):
        self.ntype = 'python'
        name = os.path.basename(fname)
        self.name = name.replace('.m', '')
        self.nb = Notebook()
        self.fname = fname
        self._excercise_num = 1
        self.excercises = []

    def convert(self, out_dir='.', ntype='python'):
        self.ntype = ntype.lower()
        with open(self.fname, 'rb') as fid:
            text = fid.read().decode('utf-8', 'replace')

        text = text.replace('\ufffd', ' ')
        lines = text.splitlines()

        header = lines[0][3:].rstrip()
        header = [header, '=' * len(header)]
        header += [lib.INTRO % ntype]
        self.nb.add_markdown(header + lib.MATH_CMDS)

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
        basename = os.path.basename(fname)
        path = os.path.join(out_dir, basename)

        self.nb.save(path)

        if self.ntype == 'python':
            self._write_exercises(out_dir)
        elif self.ntype == 'matlab':
            convert_to_matlab_kernel(path)

    def parse_line(self, line, state):
        new_line = line
        new_state = state

        if state == 'excercise':
            if new_line.startswith('%EXO'):
                new_state = 'markdown'
                new_line = ''

        elif state == 'comment':
            if new_line.startswith('%CMT'):
                new_state = 'markdown'
                new_line = ''

        elif re.match(lib.SECTION_HEADING, new_line):
            new_state = 'markdown'
            new_line = new_line[3:]
            new_line += '\n' + '-' * len(new_line)

        elif new_line.startswith('perform_toolbox_installation('):
            new_state = 'install'
            # get the package names
            new_line = re.split("'(\w*?)'", new_line)
            new_line = ' '.join(new_line[1::2])

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
            pass

        elif new_state == 'code' and self.ntype == 'python':
            new_line = self.parse_code(new_line)

        return new_state, new_line

    def parse_code(self, line):
        if line.rstrip().endswith(';'):
            line = line.rstrip()[:-1]

        for (pattern, repl) in lib.PY_REPLS:
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

        if state == 'excercise':
            header = '__Exercise %s__' % self._excercise_num
            lines = [header, '']

            code_lines = []
            for line in out_lines:
                if line.startswith('%'):
                    lines.append(line[3:])
                else:
                    code_lines.append(line)
            lines = self._parse_markdown(lines)

            self.excercises.append((lines[2:], code_lines))

            nb.add_markdown(lines)

            if self.ntype == 'python':
                nb.add_code('solutions.exo%s()' % self._excercise_num)
            else:
                nb.add_code('%%%%matlab\nexo%s()' % self._excercise_num)

            self._excercise_num += 1
            if self.ntype == 'matlab':
                nb.add_code("%%matlab\n%% Insert your code here.")
            else:
                nb.add_code("## Insert your code here.")

        elif state == 'markdown':
            nb.add_markdown(self._parse_markdown(out_lines))

        elif state == 'code':
            if not self.ntype == 'python':
                out_lines.insert(0, '%%{0}'.format(self.ntype))
            nb.add_code(out_lines)

        elif state == 'install':
            func = getattr(self, 'get_%s_intro' % self.ntype)
            func(out_lines[0].split())

    def get_python_intro(self, toolboxes):
        setup = r"""
        from __future__ import division
        import nt_toolbox as nt
        from nt_solutions import {0} as solutions
        %matplotlib inline
        %load_ext autoreload
        %autoreload 2
        """.format(self.name)

        self.nb.add_code(self._reformat(setup))

    def get_matlab_intro(self, toolboxes):
        setup = ['%load_ext pymatbridge']
        self.nb.add_code(setup)

        setup = ['%%matlab']
        for toolbox in toolboxes:
            setup += ["addpath('toolbox_%s')" % toolbox]
        setup += ["addpath('solutions/%s')" % self.name]
        self.nb.add_code(setup)

    def get_scilab_intro(self, toolboxes):
        setup = ['%load_ext scilab2py.ipython']
        self.nb.add_code(setup)

        setup = ['%%scilab']
        for toolbox in toolboxes:
            setup += ["addpath('toolbox_%s')" % toolbox]
        setup += ["addpath('solutions/%s')" % self.name]
        self.nb.add_code(setup)

    def _parse_markdown(self, lines):
        lines = self._handle_links(lines)
        return self._handle_latex(lines)

    @staticmethod
    def _handle_links(lines):
        links = []
        new_lines = []
        for line in lines:
            matches = re.findall(lib.LINK, line)
            for match in matches:
                if not isinstance(match, tuple):
                    continue
                link, title = match
                links.append(link[1:-1])
                line = line.replace(link, '')
                new_link = '[%s][%s]' % (title[1:-2], len(links))
                line = line.replace(title, new_link)
            biblio_links = re.findall(lib.BIBLIO_LINK, line)
            for link in biblio_links:
                line = line.replace('<#biblio %s>' % link,
                                    '%s(#biblio)' % link)
            new_lines.append(line)
        if links:
            new_lines.append('')
        for (ind, link) in enumerate(links):
            new_lines.append('[%s]:%s' % (ind + 1, link))

        return new_lines

    @staticmethod
    def _handle_latex(lines):
        output = []
        for line in lines:
            while line.startswith('%'):
                line = line[1:]
            if line.startswith(' '):
                line = line[1:]
            if '\\' in line:
                for (pattern, repl) in lib.MATH_REPLS:
                    line = re.sub(pattern, repl, line)
            output.append(line.rstrip())
        return output

    @staticmethod
    def _reformat(text):
        lines = text.splitlines()
        return '\n'.join([l.strip() for l in lines])

    def _write_exercises(self, out_dir):
        sfile = 'solutions/%s.py' % self.name
        sfile = os.path.join(out_dir, sfile)
        with open(sfile, 'w') as fid:
            for (ind, (comments, lines)) in enumerate(self.excercises):
                fid.write('def exo%s():\n    """\n' % (ind + 1))
                for comment in comments:
                    fid.write('    %s\n' % comment)
                fid.write('    """\n')
                for line in lines:
                    line = self.parse_code(line)
                    if line.strip():
                        fid.write('    %s\n' % line)
                fid.write('\n\n')

if __name__ == '__main__':
    usage = 'm2nb_converter.py fname [destination_dir]'
    ntype = 'python'
    ntype = 'matlab'
    if len(sys.argv) >= 3:
        Converter(sys.argv[1]).convert(sys.argv[2], ntype)
    elif len(sys.argv) >= 2:
        Converter(sys.argv[1]).convert('.', ntype)
    else:
        print(usage)
