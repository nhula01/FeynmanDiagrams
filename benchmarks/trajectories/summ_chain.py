import json
d = json.load(open('chain8_contour.json'))
for k, v in d['results'].items():
    o = v['diagrams_maxdeg5']['orders']
    ot = {q: (t['N_star'], round(t['value'], 4), round(t['delta'], 4)) for q, t in v.get('optimal_truncation_maxdeg5', {}).items()}
    print(k, 'N3', [round(x, 4) for x in o[3]], 'N4', [round(x, 4) for x in o[4]], 'N6', [round(x, 4) for x in o[6]],
          'G', [round(x, 4) for x in v['gaussian']['value']], 'MF', round(v['mean_field'][0], 4), ot)

