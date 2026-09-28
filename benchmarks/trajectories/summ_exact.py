import json, sys, glob
for f in sys.argv[1:]:
    for f2 in sorted(glob.glob(f)):
        print('##', f2)
        for l in open(f2):
            try: d = json.loads(l)
            except Exception: print('  (partial line)'); continue
            th = {k: [round(x, 5) for x in v['ratios'][:3]] for k, v in d.items() if k.startswith('theta_r')}
            fw = {k: [round(x, 5) for x in v[1:]] for k, v in d['finite_window'].items()}
            print(' cut', d['cut'], 'n', [round(x, 5) for x in d['n']], 'top %.1e' % d['top'][0], th, 'fd', [round(x, 5) for x in d['fd_h0.05']],
                  'gap', d['L0_eigs'][1], 'wind', d.get('winding_max'), 'wall %.0f' % d['wall'])
            print('   fw', fw)

