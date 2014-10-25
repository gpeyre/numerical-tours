import json
import re

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

MATH_REPLS = [(re.compile(r'\\\['), '$$'),  # replace latex delimiters
              (re.compile(r'\\\]'), '$$'),
              (re.compile(r'\\\('), '$'),
              (re.compile(r'\\\)'), '$'),
              ]


class Notebook(dict):
    """An IPython Notebook builder tailored for numerical-tours.
    """
    def __init__(self, name):
        super(Notebook, self).__init__()
        self._name = name
        self.update(dict(metadata=dict(name=""),
                         nbformat=3,
                         nbformat_minor=0,
                         worksheets=[dict(cells=[])]))
        self.add_markdown(LATEX_COMMANDS)

    def add_heading(self, source, level=1):
        source = self._handle_items(source)
        h = dict(cell_type="heading",
                 level=level,
                 metadata=dict(),
                 source=source)
        self['worksheets'][0]['cells'].append(h)

    def add_markdown(self, source):
        source = self._handle_items(source)
        output = []
        for line in source:
                while line.startswith('%'):
                    line = line[1:]
                if line.startswith(' '):
                    line = line[1:]
                if '\\' in line:
                    for (pattern, repl) in MATH_REPLS:
                        line = re.sub(pattern, repl, line)
                output.append(line.rstrip())
        md = dict(cell_type="markdown",
                  metadata={},
                  source=output)
        if output:
            self['worksheets'][0]['cells'].append(md)

    def add_code(self, source, outputs=None):
        outputs = self._handle_items(outputs)
        source = self._handle_items(source)
        code = dict(cell_type="code",
                    collapsed=False,
                    input=source,
                    language="python",
                    metadata={},
                    outputs=outputs)
        self['worksheets'][0]['cells'].append(code)

    def save(self, path):
        with open(path, 'w') as fid:
            json.dump(self, fid, indent=2, sort_keys=True)

    @staticmethod
    def _handle_items(items):
        if items is None:
            items = []
        elif not isinstance(items, list):
            items = [str(items)]
        elif items:
            new_items = []
            for item in items:
                item = str(item)
                if not item.endswith('\n'):
                    item += '\n'
                new_items.append(item)
            items = new_items
        # remove leading or trailing empty whitespace
        if items:
            items = ''.join(items).strip().splitlines()
        return items


if __name__ == '__main__':
    nb = Notebook()
    nb.add_heading("Hello, World!")
    nb.add_markdown('Introduction to Hello World!')
    # the intro markdown block was added here, with the latex commands
    nb.add_markdown(r'$x = \la * \ga$')
    # this will be replaced by the code intro
    nb.add_code("")
    nb.add_excercise('Explore the affect of $\sigma$')
    nb.save('test.ipynb')
