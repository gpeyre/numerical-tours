
fname = "../python/inverse_2_deconvolution_variational.ipynb"
lines = []
with open(fname) as fid:
    for line in fid:
        if 'newcommand' in line:
            if line.strip().startswith('"\\\\newcommand'):
                line = line.replace('"\\\\newcommand', '"$\\\\newcommand')
            if line.strip().startswith('"\\\\renewcommand'):
                line = line.replace('"\\\\renewcommand', '"$\\\\newcommand')
            line = line.replace('\\n"', '$\\n"')
            print(line)
        if not line.strip() == '"$"':
            lines.append(line)


with open(fname, 'w') as fid:
    fid.writelines(lines)
