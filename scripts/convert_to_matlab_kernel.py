import json


kernelspec = """{
 "metadata": {
  "kernelspec": {
   "display_name": "Matlab",
   "language": "matlab",
   "name": "matlab_kernel"
  },
  "language_info": {
   "file_extension": ".m",
   "help_links": [
    {
     "text": "MetaKernel Magics",
     "url": "https://github.com/calysto/metakernel/blob/master/metakernel/magics/README.md"
    }
   ],
   "mimetype": "text/x-matlab",
   "name": "matlab"
  }
}
}"""


def convert_to_matlab_kernel(fname):

    print('Converting %s' % fname)
    with open(fname) as fid:
        data = json.load(fid)

    data['metadata'] = json.loads(kernelspec)['metadata']

    cells = data['worksheets'][0]['cells']

    cells = [c for c in cells if not
             (c['cell_type'] == 'code' and c['input'] and
              'pymatbridge' in c['input'][0])]

    cells = [c for c in cells if not
             (c['cell_type'] == 'markdown'
              and c['source'][0].startswith('Installation'))]

    for cell in cells:

        if cell['cell_type'] != 'code':
            continue

        if not cell['input']:
            continue

        if cell['input'][0].startswith('%%matlab'):
            cell['input'] = cell['input'][1:]

    data['worksheets'][0]['cells'] = cells

    with open(fname, 'w') as fid:
        json.dump(data, fid)


if __name__ == '__main__':
    import sys
    import os
    fname = '../matlab/audio_1_processing.ipynb'

    if '--batch' in sys.argv:
        matlab_dir = '../matlab'
        files = os.listdir(matlab_dir)
        for fname in files:
            if fname.endswith('.ipynb'):
                path = os.path.join(matlab_dir, fname)
                convert_to_matlab_kernel(path)

    elif len(sys.argv) > 1:
        fname = sys.argv[1]
        convert_to_matlab_kernel(fname)
