import os
import re
import sys

from nb_template import Notebook

SECTION_HEADING = re.compile('\A%% \w')

PY_REPLS = [(re.compile('@\((.*?)\)'),  # replace anon funcs with lambdas
             lambda m: 'lambda %s: ' % m.groups()[0]),
            (re.compile('\A\W*?end\W*?\Z'), ''),  # kill "end" lines
            (re.compile('\A\W*?clf\W*?\Z'), ''),  # kill "clf" lines
            ]


GITHUB_LINK = 'https://github.com/gpeyre/numerical-tours/archive/master.zip'
IPYTHON_LINK = 'http://ipython.org/install.html'
INSTALLATION = """
Installation
------------
You need to download [numerical_tours](%s)
and have the IPython notebook [installed](%s) to run the code
""" % (GITHUB_LINK, IPYTHON_LINK)

MATH_REPLS = [(re.compile(r'\\\['), '$$'),  # replace latex delimiters
              (re.compile(r'\\\]'), '$$'),
              (re.compile(r'\\\('), '$'),
              (re.compile(r'\\\)'), '$'),
              ]

LINK = re.compile(r"(\<http.*? )(_.*?_\>)")
BIBLIO_LINK = re.compile(r'\<#biblio (\[.*?\])\>')


LATEX_COMMANDS = r"""
$\newcommand{\dotp}[2]{\langle #1, #2 \rangle}
\newcommand{\enscond}[2]{\lbrace #1, #2 \rbrace}
\newcommand{\pd}[2]{ \frac{ \partial #1}{\partial #2} }
\newcommand{\umin}[1]{\underset{#1}{\min}\;}
\newcommand{\umax}[1]{\underset{#1}{\max}\;}
\newcommand{\umin}[1]{\underset{#1}{\min}\;}
\newcommand{\uargmin}[1]{\underset{#1}{argmin}\;}
\newcommand{\norm}[1]{\|#1\|}
\newcommand{\abs}[1]{\left|#1\right|}
\newcommand{\choice}[1]{ \left\{  \begin{array}{l} #1 \end{array} \right. }
\newcommand{\pa}[1]{\left(#1\right)}
\newcommand{\diag}[1]{{diag}\left( #1 \right)}
\newcommand{\qandq}{\quad\text{and}\quad}
\newcommand{\qwhereq}{\quad\text{where}\quad}
\newcommand{\qifq}{ \quad \text{if} \quad }
\newcommand{\qarrq}{ \quad \Longrightarrow \quad }
\newcommand{\ZZ}{\mathbb{Z}}
\newcommand{\CC}{\mathbb{C}}
\newcommand{\RR}{\mathbb{R}}
\newcommand{\EE}{\mathbb{E}}
\newcommand{\Zz}{\mathcal{Z}}
\newcommand{\Ww}{\mathcal{W}}
\newcommand{\Vv}{\mathcal{V}}
\newcommand{\Nn}{\mathcal{N}}
\newcommand{\NN}{\mathcal{N}}
\newcommand{\Hh}{\mathcal{H}}
\newcommand{\Bb}{\mathcal{B}}
\newcommand{\Ee}{\mathcal{E}}
\newcommand{\Cc}{\mathcal{C}}
\newcommand{\Gg}{\mathcal{G}}
\newcommand{\Ss}{\mathcal{S}}
\newcommand{\Pp}{\mathcal{P}}
\newcommand{\Ff}{\mathcal{F}}
\newcommand{\Xx}{\mathcal{X}}
\newcommand{\Mm}{\mathcal{M}}
\newcommand{\Ii}{\mathcal{I}}
\newcommand{\Dd}{\mathcal{D}}
\newcommand{\Ll}{\mathcal{L}}
\newcommand{\Tt}{\mathcal{T}}
\newcommand{\si}{\sigma}
\newcommand{\al}{\alpha}
\newcommand{\la}{\lambda}
\newcommand{\ga}{\gamma}
\newcommand{\Ga}{\Gamma}
\newcommand{\La}{\Lambda}
\newcommand{\si}{\sigma}
\newcommand{\Si}{\Sigma}
\newcommand{\be}{\beta}
\newcommand{\de}{\delta}
\newcommand{\De}{\Delta}
\renewcommand{\phi}{\varphi}
\renewcommand{\th}{\theta}
\newcommand{\om}{\omega}
\newcommand{\Om}{\Omega}
$
""".strip().splitlines()


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
        self.name = name.replace('.m', '')
        self.nb = Notebook()
        self.fname = fname
        self._excercise_num = 1

    def convert(self):
        with open(self.fname, 'rb') as fid:
            text = fid.read().decode('utf-8', 'replace')

        text = text.replace('\ufffd', ' ')
        lines = text.splitlines()

        header = lines[0][3:].rstrip()
        header = [header, '=' * len(header)]
        self.nb.add_markdown(header + LATEX_COMMANDS)

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
            path = os.path.join(dname, 'notebooks', fname)

        self.nb.save(path)

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

        elif re.match(SECTION_HEADING, new_line):
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

        if state == 'excercise':
            header = 'Exercise %s' % self._excercise_num

            # TODO: write out appropriate execercise here
            # for ../python/numerical_tours/solutions/%name.py
            # ../matalb/solutions/%
            """
            if new_line.startswith('%'):
                new_state = 'excercise'
                new_line = new_line[3:]

            else:
                new_state = 'excercise'
                new_line = ''

            out_lines = [header, '-' * len(header)] + out_lines
            nb.add_markdown(self._parse_markdown(out_lines))
            """

            if self.ntype == 'python':
                nb.add_code('excercises.ex%s()' % self._excercise_num)
            else:
                nb.add_code('%%%%matlab\nex%s()' % self._excercise_num)

            self._excercise_num += 1
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

    def get_python_intro(self, *toolboxes):
        setup = r"""
        from __future__ import division
        import .nt_toolbox as nt
        from .solutions import {0} as exercises
        %matplotlib inline
        %load_ext autoreload
        %autoreload 2
        """.format(self.name)

        self.nb.add_code(self._reformat(setup))
        self.nb.add_markdown(INSTALLATION)

    def get_matlab_intro(self, toolboxes):
        setup = ['%load_ext pymatbridge']
        self.nb.add_code(setup)

        setup = ['%%matlab']
        for toolbox in toolboxes:
            setup += ["addpath('../toolbox_%s')" % toolbox]
        setup += ["addpath('../solutions/%s')" % self.name]
        self.nb.add_code(setup)

        notice = """You must also install `pymatbridge`:

        ```
        pip install pymatbridge
        ```
        """
        notice = INSTALLATION + notice
        self.nb.add_markdown(self._reformat(notice))

    def get_scilab_intro(self, toolboxes):
        setup = ['%load_ext scilab2py.ipython']
        self.nb.add_code(setup)

        setup = ['%%scilab']
        for toolbox in toolboxes:
            setup += ["addpath('../toolbox_%s')" % toolbox]
        setup += ["addpath('../solutions/%s')" % self.name]
        self.nb.add_code(setup)

        notice = """You must also install `scilab2py`:

        ```
        pip install scilab2py
        ```
        """
        notice = INSTALLATION + notice
        self.nb.add_markdown(self._reformat(notice))

    def _parse_markdown(self, lines):
        lines = self._handle_links(lines)
        return self._handle_latex(lines)

    @staticmethod
    def _handle_links(lines):
        links = []
        new_lines = []
        for line in lines:
            matches = re.findall(LINK, line)
            for match in matches:
                if not isinstance(match, tuple):
                    continue
                link, title = match
                links.append(link[1:-1])
                line = line.replace(link, '')
                new_link = '[%s][%s]' % (title[1:-2], len(links))
                line = line.replace(title, new_link)
            biblio_links = re.findall(BIBLIO_LINK, line)
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
                for (pattern, repl) in MATH_REPLS:
                    line = re.sub(pattern, repl, line)
            output.append(line.rstrip())
        return output

    @staticmethod
    def _reformat(text):
        lines = text.splitlines()
        return '\n'.join([l.strip() for l in lines])

if __name__ == '__main__':
    if len(sys.argv) >= 2:
        fname = sys.argv[1]
    else:
        fname = '../matlab/meshwav_1_subdivision_curves.m'
    Converter(fname, 'python').convert()
