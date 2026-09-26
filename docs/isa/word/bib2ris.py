"""refs.bib -> references_zotero.ris (RIS for Zotero, UTF-8, CRLF line endings).

usage (from docs/isa):  python word/bib2ris.py
Type map: article JOUR, inproceedings CPAPER, book BOOK, misc (data sheets) ELEC.
LaTeX accents and math are converted to Unicode; the BibTeX key is kept in ID.
"""
from __future__ import annotations

import os
import re

from pylatexenc.latex2text import LatexNodes2Text

HERE = os.path.dirname(os.path.abspath(__file__))
ISA = os.path.dirname(HERE)
L2T = LatexNodes2Text(math_mode='text')
TYPES = {'article': 'JOUR', 'inproceedings': 'CPAPER', 'book': 'BOOK', 'misc': 'ELEC'}


def entries(text):
    for m in re.finditer(r'@(\w+)\{([^,]+),', text):
        i, d = m.end(), 1
        j = i
        while d:
            d += {'{': 1, '}': -1}.get(text[j], 0)
            j += 1
        body = text[i:j - 1]
        fields, k = {}, 0
        for fm in re.finditer(r'(\w+)\s*=\s*\{', body):
            if fm.start() < k:
                continue
            a, dd = fm.end(), 1
            b = a
            while dd:
                dd += {'{': 1, '}': -1}.get(body[b], 0)
                b += 1
            fields[fm.group(1).lower()] = body[a:b - 1]
            k = b
        yield m.group(1).lower(), m.group(2).strip(), fields


def txt(s):
    s = s.replace('$H_\\infty$', 'H∞').replace('{$H_\\infty$}', 'H∞')
    out = L2T.latex_to_text(s)
    return re.sub(r'\s+', ' ', out).strip()


def main():
    bib = open(os.path.join(ISA, 'refs.bib'), encoding='utf-8').read()
    recs = []
    for typ, key, f in entries(bib):
        r = [('TY', TYPES[typ])]
        for a in re.split(r'\s+and\s+', f.get('author', '')):
            a = a.strip()
            if not a:
                continue
            if a.startswith('{') and a.endswith('}') and a.count('{') == 1:
                r.append(('AU', txt(a) + ','))             # corporate author: kept whole by Zotero
            else:
                r.append(('AU', txt(a)))
        r.append(('TI', txt(f['title'])))
        if 'journal' in f:
            r += [('T2', txt(f['journal'])), ('JO', txt(f['journal']))]
        if 'booktitle' in f:
            r.append(('T2', txt(f['booktitle'])))
        if 'year' in f:
            r += [('PY', f['year']), ('DA', f['year'])]
        for bk, rk in (('volume', 'VL'), ('number', 'IS'), ('publisher', 'PB'), ('address', 'CY'),
                       ('edition', 'ET')):
            if bk in f:
                r.append((rk, txt(f[bk])))
        if 'pages' in f:
            p = f['pages'].replace('--', '-').split('-')
            r.append(('SP', p[0]))
            if len(p) > 1:
                r.append(('EP', p[1]))
        if 'doi' in f:
            r += [('DO', f['doi']), ('UR', 'https://doi.org/' + f['doi'])]
        elif 'url' in f:
            r.append(('UR', f['url']))
        if 'note' in f:
            r.append(('N1', txt(f['note'])))
        if typ == 'misc':
            r.append(('M3', 'Data sheet'))
        r += [('ID', key), ('LA', 'en'), ('ER', '')]
        recs.append(r)
    out = os.path.join(ISA, 'references_zotero.ris')
    with open(out, 'w', encoding='utf-8', newline='\r\n') as fh:
        for r in recs:
            for k, v in r:
                fh.write(f'{k}  - {v}\n')
            fh.write('\n')
    print(len(recs), 'references ->', out)


if __name__ == '__main__':
    main()
