import os
import re
import sys
try:
    from smop.main import main as smop_main
except ImportError:
    print('Please `pip install smop`')


def compile(fname):
    """Compile an M-file to a Python file using smop
    """
    with open(fname) as fid:
        data = fid.read()
    temp_mfile = fname.replace('.m', '_temp.m')
    with open(temp_mfile, 'wb') as fid:
        fid.write('function foo()\n')
        fid.write(data)
        fid.write('\nend')
    pyfile = fname.replace('.m', '.py')
    sys.argv = ['', temp_mfile, '-o', pyfile, '-v']
    smop_main()
    with open(pyfile, 'rb') as fid:
        data = fid.readlines()
    data = data[8:]
    output = []
    for line in data:
        line = line[4:]
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
        output.append(line)
    with open(pyfile, 'wb') as fid:
        fid.write(''.join(output))
    os.remove(temp_mfile)


if __name__ == '__main__':
    if not len(sys.argv) >= 2:
        fname = 'audio_1_processing.m'
    else:
        fname = sys.argv[1]
    compile(fname)
