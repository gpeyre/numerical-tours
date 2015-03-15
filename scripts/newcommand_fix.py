import glob
import json

for fname in glob.glob('../python/todo/*.ipynb'):
    print(fname)
    lines = []
    with open(fname) as fid:
        try:
            tree = json.load(fid)
        except (ValueError, TypeError) as e:
            print(e)
            continue

    for cell in tree['worksheets'][0]['cells']:
        lines = []
        if 'source' in cell:
            for line in cell['source']:
                if line.startswith('\\newcommand'):
                    line = '$' + line
                if line.startswith('\\renewcommand'):
                    line = '$' + line.replace('\\renewcommand', '\\newcommand')
                line = line.replace('}\n', '}$\n')
                if 'newcommand' in line:
                    if not line.count('$') == 2:
                        line += '$'
                if not line == '$':
                    lines.append(line)
            cell['source'] = lines

    with open(fname, 'w') as fid:
        print('saving', fname)
        json.dump(tree, fid)
