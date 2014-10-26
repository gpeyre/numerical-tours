import json


class Notebook(dict):

    """An IPython Notebook builder tailored for numerical-tours.
    """

    def __init__(self):
        super(Notebook, self).__init__()
        self.update(dict(metadata=dict(name=""),
                         nbformat=3,
                         nbformat_minor=0,
                         worksheets=[dict(cells=[])]))

    def add_markdown(self, source):
        source = self._handle_items(source)
        md = dict(cell_type="markdown",
                  metadata={},
                  source=source)
        if source:
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
    nb.add_markdown(r'$x = \la * \ga$')
    nb.add_code('a = 1')
    nb.save('test.ipynb')
