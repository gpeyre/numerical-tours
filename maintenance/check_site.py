"""Validate the generated public catalogue and every exported Python tour."""
from html.parser import HTMLParser
import json
from pathlib import Path
import sys
from urllib.parse import urlparse, unquote

class Elements(HTMLParser):
    def __init__(self):
        super().__init__()
        self.tags = []
    def handle_starttag(self, tag, attrs):
        self.tags.append((tag, dict(attrs)))

root = Path(sys.argv[1]).resolve()
source = Path(__file__).resolve().parents[1]
tours = json.loads((source / '_data/python_catalog.json').read_text())
assert len(tours) == 58
for route in ['index.html', 'python/index.html', 'archive/index.html', 'installation_python/index.html', 'about/index.html']:
    path = root / route
    assert path.is_file(), route
    parsed = Elements(); parsed.feed(path.read_text())
    assert any(a.get('id') == 'main' for _, a in parsed.tags), f'Missing main landmark: {route}'
    for tag, attrs in parsed.tags:
        for key in ['href', 'src']:
            url = attrs.get(key, '')
            if not url.startswith('/') or url.startswith('//'):
                continue
            target = root / unquote(urlparse(url).path).lstrip('/')
            if target.is_dir():target = target / 'index.html'
            assert target.exists(), (route, url)
for tour in tours:
    page = root / tour['url'].strip('/') / 'index.html'
    assert page.is_file()
    assert (root / tour['image'].lstrip('/')).stat().st_size > 1000
    assert (root / 'downloads' / (tour['slug']+'.ipynb')).is_file()
    text = page.read_text()
    parsed = Elements(); parsed.feed(text)
    headings = [(tag, attrs['id']) for tag, attrs in parsed.tags if tag in ['h1', 'h2', 'h3']]
    assert sum(tag == 'h1' for tag, _ in headings) == 1, tour['slug']
    heading_ids = [anchor for _, anchor in headings]
    assert len(heading_ids) == len(set(heading_ids)), tour['slug']
    contents = Elements(); contents.feed(text.split('<aside class="notebook-sidebar"', 1)[1].split('</aside>', 1)[0])
    targets = [unquote(attrs['href'][1:]) for tag, attrs in contents.tags if tag == 'a']
    assert targets == heading_ids, (tour['slug'], 'Contents must link to every heading in order')
    assert '/assets/css/notebook.css' in text and '/assets/js/notebook.js' in text
    assert 'Open in Colab' in text and tour['colab'] in text
    assert 'References' in text
    assert 'http://cdn.mathjax.org/' not in text
    notebook = json.loads((root / 'downloads' / (tour['slug']+'.ipynb')).read_text())
    assert not any(o.get('output_type') == 'error' for c in notebook['cells'] for o in c.get('outputs', []))
    assert all(c['execution_count'] is not None for c in notebook['cells'] if c['cell_type'] == 'code')
print(f'Passed: {len(tours)} reading pages, section navigation, previews, downloads, Colab links, and primary navigation routes.')
