import glob
import json

for fname in glob.glob('../matlab/*.ipynb'):
    print(fname)
    lines = []
    with open(fname) as fid:
        
        """
        for line in fid:
            if 'newcommand' in line: 
                from PyQt4.QtCore import pyqtRemoveInputHook; pyqtRemoveInputHook()
                import ipdb; ipdb.set_trace()
                pass
                
                                  
                line = line.replace('\$\newcommand{', '$\\newcommand{')
                #line = line.replace('}$"', '}$",')
            if not line.strip() == '"$"':
                lines.append(line)
        """
        try:
            tree = json.load(fid)
        except (ValueError, TypeError) as e:
            print(e)
            from PyQt4.QtCore import pyqtRemoveInputHook; pyqtRemoveInputHook()
            import ipdb; ipdb.set_trace()
            pass
            
            continue

    for line in tree['worksheets'][0]['cells'][0]['source']:
        if line.startswith('\\newcommand'):
            line = '$' + line
        if line.startswith('\\renewcommand'):
            line = '$' + line.replace('\\renewcommand', '\\newcommand')
        line = line.replace('}\n', '}$\n')
            
        if not line == '$':
            lines.append(line)

    tree['worksheets'][0]['cells'][0]['source'] = lines
    
    with open(fname, 'w') as fid:
        #fid.writelines(lines)
        print('saving', fname)
        json.dump(tree, fid)
