import re


PY_REPLS = [(re.compile('@\((.*?)\)'),  # replace anon funcs with lambdas
             lambda m: 'lambda %s: ' % m.groups()[0]),
            (re.compile('\A\W*?end\W*?\Z'), ''),  # kill "end" lines
            (re.compile('\A\W*?clf\W*?\Z'), ''),  # kill "clf" lines
            (re.compile('(exo\d+)'),
             lambda m: 'solutions.%s' % m.groups()[0])
            ]

MAT_REPLS = [(re.compile('lambda (.*?):'),  # replace lambda funcs with anons
              lambda m: '@(%s) ' % m.groups()[0]),
             (re.compile('solutions.(exo\d+)'),
              lambda m: '%s' % m.groups()[0])
             ]

MATH_REPLS = [(re.compile(r'\\\['), '$$'),  # replace latex delimiters
              (re.compile(r'\\\]'), '$$'),
              (re.compile(r'\\\('), '$'),
              (re.compile(r'\\\)'), '$'),
              ]

LINK = re.compile(r"(\<http.*? )(_.*?_\>)")
BIBLIO_LINK = re.compile(r'\<#biblio (\[.*?\])\>')
SECTION_HEADING = re.compile('\A%% \w')

INTRO = """
*Important:* Please read the [installation page](http://gpeyre.github.io/numerical-tours/installation_%s/) for details about how to install the toolboxes.
"""

MATH_CMDS = r"""
{\dotp}[2]{\langle #1, #2 \rangle}
{\enscond}[2]{\lbrace #1, #2 \rbrace}
{\pd}[2]{ \frac{ \partial #1}{\partial #2} }
{\umin}[1]{\underset{#1}{\min}\;}
{\umax}[1]{\underset{#1}{\max}\;}
{\umin}[1]{\underset{#1}{\min}\;}
{\uargmin}[1]{\underset{#1}{argmin}\;}
{\norm}[1]{\|#1\|}
{\abs}[1]{\left|#1\right|}
{\choice}[1]{ \left\{  \begin{array}{l} #1 \end{array} \right. }
{\pa}[1]{\left(#1\right)}
{\diag}[1]{{diag}\left( #1 \right)}
{\qandq}{\quad\text{and}\quad}
{\qwhereq}{\quad\text{where}\quad}
{\qifq}{ \quad \text{if} \quad }
{\qarrq}{ \quad \Longrightarrow \quad }
{\ZZ}{\mathbb{Z}}
{\CC}{\mathbb{C}}
{\RR}{\mathbb{R}}
{\EE}{\mathbb{E}}
{\Zz}{\mathcal{Z}}
{\Ww}{\mathcal{W}}
{\Vv}{\mathcal{V}}
{\Nn}{\mathcal{N}}
{\NN}{\mathcal{N}}
{\Hh}{\mathcal{H}}
{\Bb}{\mathcal{B}}
{\Ee}{\mathcal{E}}
{\Cc}{\mathcal{C}}
{\Gg}{\mathcal{G}}
{\Ss}{\mathcal{S}}
{\Pp}{\mathcal{P}}
{\Ff}{\mathcal{F}}
{\Xx}{\mathcal{X}}
{\Mm}{\mathcal{M}}
{\Ii}{\mathcal{I}}
{\Dd}{\mathcal{D}}
{\Ll}{\mathcal{L}}
{\Tt}{\mathcal{T}}
{\si}{\sigma}
{\al}{\alpha}
{\la}{\lambda}
{\ga}{\gamma}
{\Ga}{\Gamma}
{\La}{\Lambda}
{\si}{\sigma}
{\Si}{\Sigma}
{\be}{\beta}
{\de}{\delta}
{\De}{\Delta}
{\phi}{\varphi}
{\th}{\theta}
{\om}{\omega}
{\Om}{\Omega}
{\eqdef}{\equiv}
""".strip().splitlines()

MATH_CMDS = '$\\newcommand' + '$\n$\\newcommand'.join(MATH_CMDS) + '$'
MATH_CMDS = MATH_CMDS.splitlines()
