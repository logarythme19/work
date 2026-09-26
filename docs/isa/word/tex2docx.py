"""LaTeX manuscript -> Word (.docx) with native equations, fields and bookmarks.

usage (from docs/isa):  python word/tex2docx.py
Needs: pandoc >= 3, a compiled main.aux (numbers are taken from LaTeX itself),
word/png/*.png (figures), word/elsevier-harvard.csl.

Pipeline
 1. flatten main.tex (\\input), take every label number from main.aux;
 2. rewrite what pandoc cannot number: display equations are followed by an
    [[EQ]] marker, captions start with a [[CAP]] marker, \\ref/\\eqref become
    [[REF]] / [[LNK]] markers, \\resizebox wrappers and run-in \\paragraph
    headings are unwrapped, siunitx inside math is expanded;
 3. pandoc -> docx (OMML equations, citeproc author-year with links);
 4. post-process word/document.xml: each equation becomes a borderless
    two-cell table (equation centred, "(n)" right, SEQ Equation field inside a
    bookmark); captions get "Figure/Table" + SEQ field inside a bookmark;
    cross-references become REF fields (\\h) with the LaTeX number as cached
    result; section references become internal hyperlinks; a PAGE footer is
    added. Word shows the cached numbers and renumbers everything on F9.
"""
from __future__ import annotations

import copy
import os
import re
import subprocess
import sys
import zipfile

from lxml import etree

HERE = os.path.dirname(os.path.abspath(__file__))
ISA = os.path.dirname(HERE)
OUT = os.path.join(HERE, 'CMO_LQI_PIDF_ISA.docx')

W = 'http://schemas.openxmlformats.org/wordprocessingml/2006/main'
M = 'http://schemas.openxmlformats.org/officeDocument/2006/math'
R = 'http://schemas.openxmlformats.org/officeDocument/2006/relationships'
NS = {'w': W, 'm': M, 'r': R}


def q(tag):
    p, t = tag.split(':')
    return '{%s}%s' % (NS[p], t)


# ---------------------------------------------------------------- flatten
def read(path):
    return open(path, encoding='utf-8').read()


def flatten(text, base):
    def rep(m):
        name = m.group(1)
        f = os.path.join(base, name if name.endswith('.tex') else name + '.tex')
        if name.startswith('figures/fig_') and os.path.exists(os.path.join(HERE, 'png', os.path.basename(name) + '.png')):
            return '\\includegraphics[width=\\textwidth]{%s}' % os.path.basename(name)
        return flatten(read(f), base)
    text = re.sub(r'(?m)(?<!\\)%.*$', '', text)          # comments
    return re.sub(r'\\input\{([^}]*)\}', rep, text)


def labels():
    aux = read(os.path.join(ISA, 'main.aux'))
    lab = dict(re.findall(r'\\newlabel\{([^}]*)\}\{\{([^}]*)\}', aux))
    return {k: v.replace('Appendix ~', '').strip() for k, v in lab.items()}


def bm(label):
    return re.sub(r'[^A-Za-z0-9_]', '_', label)[:38]


def group(s, i):
    """Return (content, end) of the brace group starting at s[i] == '{'."""
    assert s[i] == '{', s[i:i + 20]
    d = 0
    for j in range(i, len(s)):
        if s[j] == '{':
            d += 1
        elif s[j] == '}':
            d -= 1
            if d == 0:
                return s[i + 1:j], j + 1
    raise ValueError('unbalanced')


def unwrap_resizebox(s):
    while True:
        i = s.find('\\resizebox{')
        if i < 0:
            return s
        _, j = group(s, i + len('\\resizebox'))
        _, k = group(s, j)
        body, e = group(s, k)
        s = s[:i] + body + s[e:]


UNITS = {r'\ohm': r'\Omega', r'\micro': r'\mu ', r'\percent': r'\%', r'\degree': r'^{\circ}',
         r'\celsius': r'^{\circ}C'}


def unit_math(u):
    for a, b in UNITS.items():
        u = u.replace(a, b)
    return r'\mathrm{' + u.replace('.', r'\,') + '}'


def si_in_math(m):
    def si(mm):
        return unit_math(mm.group(1))

    def SI(mm):
        v = mm.group(1)
        v = re.sub(r'([0-9.]+)e([+-]?\d+)', r'\1\\times10^{\2}', v)
        u = mm.group(2)
        return v + (r'\,' + unit_math(u) if u else '')
    m = re.sub(r'\\SI\{([^}]*)\}\{([^}]*)\}', SI, m)
    return re.sub(r'\\si\{([^}]*)\}', si, m)


def math_clean(m):
    """Glyph-robust spacing and transpose for Word/LibreOffice fonts."""
    m = si_in_math(m)
    m = m.replace('\\qquad', '\\ \\ \\ ').replace('\\quad', '\\ \\ ').replace('\\;', '\\ ')
    return re.sub(r'\^\\top\b', lambda _: '^{\\mathrm{T}}', m)


def convert_tex():
    lab = labels()
    src = flatten(read(os.path.join(ISA, 'main.tex')), ISA)
    pre, body = src.split('\\begin{document}')
    body = body.split('\\end{document}')[0]
    macros = '\n'.join(l for l in read(os.path.join(ISA, 'numbers.tex')).splitlines() if l.startswith('\\newcommand'))
    # --- front matter
    title = re.search(r'\\title\{(.*?)\}\n', body, re.S).group(1)
    abstract = re.search(r'\\begin\{abstract\}(.*?)\\end\{abstract\}', body, re.S).group(1)
    kw = re.search(r'\\begin\{keyword\}(.*?)\\end\{keyword\}', body, re.S).group(1)
    kw = '; '.join(k.strip() for k in kw.split('\\sep'))
    body = re.sub(r'\\begin\{frontmatter\}.*?\\end\{frontmatter\}', '', body, flags=re.S)
    front = ('\\begin{center}\\textbf{%s}\\end{center}\n\n\\begin{center}Anonymous author(s) --- '
             'affiliation withheld for review\\end{center}\n\n\\section*{Abstract}\n%s\n\n'
             '\\textbf{Keywords:} %s\n\n' % (title, abstract.strip(), kw))
    body = front + body
    for cmd in ('\\linenumbers', '\\bibliographystyle{elsarticle-harv}', '\\bibliography{refs}'):
        body = body.replace(cmd, '')
    body = unwrap_resizebox(body)
    body = body.replace('\\footnotesize', '').replace('\\centering', '')
    # --- appendix headings
    a, b = body.split('\\appendix')
    letters = iter('ABCDEF')
    b = re.sub(r'\\section\{([^}]*)\}', lambda m: '\\section*{Appendix %s. %s}' % (next(letters), m.group(1)), b)
    body = a + b
    # --- run-in paragraph headings
    body = re.sub(r'\\paragraph\{([^}]*)\}\s*', lambda m: '\\textbf{%s.} ' % m.group(1).rstrip('.'), body)
    # --- floats: caption markers
    def float_env(m):
        env, content = m.group(1), m.group(2)
        kind = 'Figure' if env.startswith('figure') else 'Table'
        lm = re.search(r'\\label\{([^}]*)\}', content)
        if not lm:
            raise ValueError('unlabelled float: ' + content[:80])
        num = lab[lm.group(1)]
        content = content.replace(lm.group(0), '')
        i = content.find('\\caption{')
        cap, e = group(content, i + len('\\caption'))
        content = content[:i] + '\\caption{[[CAP|%s|%s|%s]] %s}' % (kind, num, bm(lm.group(1)), cap.strip()) + content[e:]
        # pandoc places the caption from the environment; keep order
        env = env.replace('table*', 'table')      # pandoc drops captions of table*
        return '\\begin{%s}%s\\end{%s}' % (env, content, env)
    body = re.sub(r'\\begin\{(figure\*?|table\*?)\}(?:\[[^\]]*\])?(.*?)\\end\{\1\}', float_env, body, flags=re.S)
    body = body.replace('\\includegraphics[width=\\columnwidth]{figures/', '\\includegraphics[width=0.75\\textwidth]{figures/')
    body = re.sub(r'\\includegraphics\[([^\]]*)\]\{(?:figures/)?([^}.]*)(?:\.pdf)?\}',
                  lambda m: '\\includegraphics[%s]{%s/png/%s.png}' % (m.group(1), HERE, m.group(2)), body)
    # --- equations
    n = [0]
    def eq(m):
        n[0] += 1
        content = m.group(1)
        lm = re.search(r'\\label\{([^}]*)\}', content)
        if lm:
            assert lab[lm.group(1)] == str(n[0]), (lm.group(1), lab[lm.group(1)], n[0])
            content = content.replace(lm.group(0), '')
            b_ = bm(lm.group(1))
        else:
            b_ = 'eq_%d' % n[0]
        content = math_clean(content.strip())
        return '\n\n\\[\n%s\n\\]\n\n[[EQ|%d|%s]]\n\n' % (content, n[0], b_)
    body = re.sub(r'\\begin\{equation\}(.*?)\\end\{equation\}', eq, body, flags=re.S)
    # siunitx inside inline math
    body = re.sub(r'\$([^$]+)\$', lambda m: '$' + math_clean(m.group(1)) + '$', body)
    # --- cross references
    def ref(m, paren):
        l = m.group(1)
        num = lab[l]
        if l.startswith(('eq:', 'fig:', 'tab:')):
            return '[[REF|%s|%s]]' % (bm(l), '(%s)' % num if paren else num)
        if l.startswith(('sec:', 'app:')):
            return '[[LNK|%s|%s]]' % (l, num)
        return num
    body = re.sub(r'\\eqref\{([^}]*)\}', lambda m: ref(m, True), body)
    body = re.sub(r'\\ref\{([^}]*)\}', lambda m: ref(m, False), body)
    body = body.replace('Appendix~[[LNK', 'Appendix [[LNK')
    doc = ('\\documentclass{article}\n\\usepackage{amsmath,siunitx}\n\\newtheorem{proposition}{Proposition}\n'
           + macros + '\n\\begin{document}\n' + body + '\n\\end{document}\n')
    path = os.path.join(HERE, 'flat.tex')
    open(path, 'w', encoding='utf-8').write(doc)
    return path, n[0]


# ------------------------------------------------------------ reference.docx
def reference_doc():
    ref = os.path.join(HERE, 'reference.docx')
    tmp = os.path.join(HERE, 'default_reference.docx')
    subprocess.run(['pandoc', '-o', tmp, '--print-default-data-file', 'reference.docx'], check=True)
    data = open(tmp, 'rb').read()
    os.remove(tmp)
    zin = zipfile.ZipFile(__import__('io').BytesIO(data))
    zout = zipfile.ZipFile(ref, 'w', zipfile.ZIP_DEFLATED)
    for it in zin.infolist():
        b = zin.read(it.filename)
        if it.filename == 'word/styles.xml':
            b = style_xml(b)
        zout.writestr(it, b)
    zout.close()
    return ref


PPR_ORDER = ['pStyle', 'keepNext', 'keepLines', 'pageBreakBefore', 'framePr', 'widowControl', 'numPr',
             'suppressLineNumbers', 'pBdr', 'shd', 'tabs', 'suppressAutoHyphens', 'kinsoku', 'wordWrap',
             'overflowPunct', 'topLinePunct', 'autoSpaceDE', 'autoSpaceDN', 'bidi', 'adjustRightInd', 'snapToGrid',
             'spacing', 'ind', 'contextualSpacing', 'mirrorIndents', 'suppressOverlap', 'jc', 'textDirection',
             'textAlignment', 'textboxTightWrap', 'outlineLvl', 'divId', 'cnfStyle', 'rPr', 'sectPr', 'pPrChange']
RPR_ORDER = ['rStyle', 'rFonts', 'b', 'bCs', 'i', 'iCs', 'caps', 'smallCaps', 'strike', 'dstrike', 'outline', 'shadow',
             'emboss', 'imprint', 'noProof', 'snapToGrid', 'vanish', 'webHidden', 'color', 'spacing', 'w', 'kern',
             'position', 'sz', 'szCs', 'highlight', 'u', 'effect', 'bdr', 'shd', 'fitText', 'vertAlign', 'rtl', 'cs',
             'em', 'lang', 'eastAsianLayout', 'specVanish', 'oMath']
STYLE_ORDER = ['name', 'aliases', 'basedOn', 'next', 'link', 'autoRedefine', 'hidden', 'uiPriority', 'semiHidden',
               'unhideWhenUsed', 'qFormat', 'locked', 'personal', 'personalCompose', 'personalReply', 'rsid', 'pPr',
               'rPr', 'tblPr', 'trPr', 'tcPr', 'tblStylePr']


def ordered_child(parent, tag, order):
    """Get or create w:<tag> inside parent at its schema position."""
    el = parent.find(q('w:' + tag))
    if el is not None:
        return el
    el = etree.Element(q('w:' + tag))
    rank = order.index(tag)
    for i, ch in enumerate(parent):
        name = etree.QName(ch).localname
        if name in order and order.index(name) > rank:
            parent.insert(i, el)
            return el
    parent.append(el)
    return el


def style_xml(b):
    t = etree.fromstring(b)

    def setp(parent, order, tag, **attrs):
        el = ordered_child(parent, tag, order)
        for k, v in attrs.items():
            el.set(q('w:' + k), v)
        return el
    dd = t.find('.//' + q('w:docDefaults'))
    rpr = dd.find('.//' + q('w:rPr'))
    setp(rpr, RPR_ORDER, 'rFonts', ascii='Times New Roman', hAnsi='Times New Roman', cs='Times New Roman',
         eastAsia='Times New Roman')
    setp(rpr, RPR_ORDER, 'sz', val='24'); setp(rpr, RPR_ORDER, 'szCs', val='24')
    for st in t.findall(q('w:style')):
        if st.get(q('w:type')) not in ('paragraph', 'character'):
            continue
        sid = st.get(q('w:styleId'))
        rp = ordered_child(st, 'rPr', STYLE_ORDER)
        for el in rp.findall(q('w:rFonts')):
            rp.remove(el)
        if sid in ('BodyText', 'FirstParagraph'):
            pp = ordered_child(st, 'pPr', STYLE_ORDER)
            setp(pp, PPR_ORDER, 'spacing', before='0', after='120', line='360', lineRule='auto')
            setp(pp, PPR_ORDER, 'jc', val='both')
        if sid.startswith('Heading') or sid == 'Title':
            setp(rp, RPR_ORDER, 'b')
            setp(rp, RPR_ORDER, 'color', val='000000')
            size = {'Heading1': '28', 'Heading2': '24', 'Heading3': '24', 'Heading4': '24'}.get(sid)
            if size:
                setp(rp, RPR_ORDER, 'sz', val=size); setp(rp, RPR_ORDER, 'szCs', val=size)
            if sid == 'Heading3':
                setp(rp, RPR_ORDER, 'i')
        if sid in ('ImageCaption', 'TableCaption', 'Caption'):
            for el in rp.findall(q('w:i')):
                rp.remove(el)
            setp(rp, RPR_ORDER, 'sz', val='20'); setp(rp, RPR_ORDER, 'szCs', val='20')
            pp = ordered_child(st, 'pPr', STYLE_ORDER)
            setp(pp, PPR_ORDER, 'jc', val='both')
        if sid == 'Compact':
            pp = ordered_child(st, 'pPr', STYLE_ORDER)
            setp(pp, PPR_ORDER, 'spacing', before='20', after='20', line='240', lineRule='auto')
            setp(pp, PPR_ORDER, 'jc', val='left')
            setp(rp, RPR_ORDER, 'sz', val='18'); setp(rp, RPR_ORDER, 'szCs', val='18')
        if sid == 'Bibliography':
            pp = ordered_child(st, 'pPr', STYLE_ORDER)
            setp(pp, PPR_ORDER, 'spacing', after='60')
            setp(pp, PPR_ORDER, 'ind', left='360', hanging='360')
        if len(rp) == 0:
            st.remove(rp)
    return etree.tostring(t, xml_declaration=True, encoding='UTF-8', standalone=True)


# ------------------------------------------------------------ post-processing
def run(text, rpr=None):
    r = etree.Element(q('w:r'))
    if rpr is not None:
        r.append(copy.deepcopy(rpr))
    t = etree.SubElement(r, q('w:t'))
    t.text = text
    t.set('{http://www.w3.org/XML/1998/namespace}space', 'preserve')
    return r


def field(instr, cached, rpr=None):
    out = []
    def fc(kind):
        r = etree.Element(q('w:r'))
        if rpr is not None:
            r.append(copy.deepcopy(rpr))
        f = etree.SubElement(r, q('w:fldChar')); f.set(q('w:fldCharType'), kind)
        return r
    out.append(fc('begin'))
    r = etree.Element(q('w:r'))
    if rpr is not None:
        r.append(copy.deepcopy(rpr))
    it = etree.SubElement(r, q('w:instrText')); it.text = ' %s ' % instr
    it.set('{http://www.w3.org/XML/1998/namespace}space', 'preserve')
    out.append(r)
    out.append(fc('separate'))
    out.append(run(cached, rpr))
    out.append(fc('end'))
    return out


class Bookmarks:
    def __init__(self, start=1000):
        self.i = start

    def wrap(self, name, elems):
        self.i += 1
        s = etree.Element(q('w:bookmarkStart')); s.set(q('w:id'), str(self.i)); s.set(q('w:name'), name)
        e = etree.Element(q('w:bookmarkEnd')); e.set(q('w:id'), str(self.i))
        return [s] + list(elems) + [e]


MARK = re.compile(r'\[\[(REF|LNK|CAP|EQ)\|([^\]]*)\]\]')


def para_text(p):
    return ''.join(t.text or '' for t in p.iter(q('w:t')))


def split_markers(p, bmk, stats):
    """Replace [[REF]], [[LNK]] and [[CAP]] markers inside the runs of paragraph p."""
    for r in list(p.iter(q('w:r'))):
        t = r.find(q('w:t'))
        if t is None or t.text is None or '[[' not in t.text:
            continue
        rpr = r.find(q('w:rPr'))
        parts, pos, new = t.text, 0, []
        for m in MARK.finditer(parts):
            if m.start() > pos:
                new.append(run(parts[pos:m.start()], rpr))
            kind, args = m.group(1), m.group(2).split('|')
            if kind == 'REF':
                new += field('REF %s \\h' % args[0], args[1], rpr); stats['ref'] += 1
            elif kind == 'LNK':
                h = etree.Element(q('w:hyperlink')); h.set(q('w:anchor'), args[0])
                h.append(run(args[1], rpr)); new.append(h); stats['lnk'] += 1
            elif kind == 'CAP':
                kindname, num, name = args
                b = etree.Element(q('w:rPr')); etree.SubElement(b, q('w:b'))
                new.append(run(kindname + ' ', b))
                if '.' in num:     # appendix numbering B.1: plain text
                    inner = [run(num, b)]
                else:
                    inner = field('SEQ %s \\* ARABIC' % kindname, num, b)
                new += bmk.wrap(name, inner)
                new.append(run(':', b))
                stats['cap'] += 1
            pos = m.end()
        if pos < len(parts):
            new.append(run(parts[pos:], rpr))
        parent = r.getparent()
        idx = parent.index(r)
        parent.remove(r)
        for k, el in enumerate(new):
            parent.insert(idx + k, el)


def eq_table(math_p, num, name, bmk, textw):
    tbl = etree.Element(q('w:tbl'))
    tp = etree.SubElement(tbl, q('w:tblPr'))
    tw = etree.SubElement(tp, q('w:tblW')); tw.set(q('w:w'), str(textw)); tw.set(q('w:type'), 'dxa')
    jc = etree.SubElement(tp, q('w:jc')); jc.set(q('w:val'), 'center')
    bd = etree.SubElement(tp, q('w:tblBorders'))
    for side in ('top', 'left', 'bottom', 'right', 'insideH', 'insideV'):
        e = etree.SubElement(bd, q('w:' + side)); e.set(q('w:val'), 'nil')
    lay = etree.SubElement(tp, q('w:tblLayout')); lay.set(q('w:type'), 'fixed')
    lk = etree.SubElement(tp, q('w:tblLook')); lk.set(q('w:val'), '0000')
    grid = etree.SubElement(tbl, q('w:tblGrid'))
    wnum = 900
    for w_ in (wnum, textw - 2 * wnum, wnum):
        g = etree.SubElement(grid, q('w:gridCol')); g.set(q('w:w'), str(w_))
    tr = etree.SubElement(tbl, q('w:tr'))
    def cell(width, content, align):
        tc = etree.SubElement(tr, q('w:tc'))
        pr = etree.SubElement(tc, q('w:tcPr'))
        cw = etree.SubElement(pr, q('w:tcW')); cw.set(q('w:w'), str(width)); cw.set(q('w:type'), 'dxa')
        va = etree.SubElement(pr, q('w:vAlign')); va.set(q('w:val'), 'center')
        p = content if content is not None else etree.Element(q('w:p'))
        ppr = p.find(q('w:pPr'))
        if ppr is None:
            ppr = etree.Element(q('w:pPr')); p.insert(0, ppr)
        for old in ppr.findall(q('w:jc')):
            ppr.remove(old)
        sp = etree.SubElement(ppr, q('w:spacing')); sp.set(q('w:before'), '60'); sp.set(q('w:after'), '60')
        j = etree.SubElement(ppr, q('w:jc')); j.set(q('w:val'), align)
        tc.append(p)
        return p
    cell(wnum, None, 'left')
    cell(textw - 2 * wnum, math_p, 'center')
    pn = cell(wnum, None, 'right')
    for el in bmk.wrap(name, [run('(')] + field('SEQ Equation \\* ARABIC', str(num)) + [run(')')]):
        pn.append(el)
    return tbl


def add_footer(z_in, files):
    """Page-number footer: footer1.xml + relationship + sectPr reference."""
    footer = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
              '<w:ftr xmlns:w="%s"><w:p><w:pPr><w:jc w:val="center"/></w:pPr>'
              '<w:r><w:fldChar w:fldCharType="begin"/></w:r><w:r><w:instrText xml:space="preserve"> PAGE </w:instrText></w:r>'
              '<w:r><w:fldChar w:fldCharType="separate"/></w:r><w:r><w:t>1</w:t></w:r><w:r><w:fldChar w:fldCharType="end"/></w:r>'
              '</w:p></w:ftr>' % W)
    files['word/footer1.xml'] = footer.encode()
    rels = files['word/_rels/document.xml.rels'].decode()
    rels = rels.replace('</Relationships>', '<Relationship Id="rIdFooter1" '
                        'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/footer" '
                        'Target="footer1.xml"/></Relationships>')
    files['word/_rels/document.xml.rels'] = rels.encode()
    ct = files['[Content_Types].xml'].decode()
    if 'footer1.xml' not in ct:
        ct = ct.replace('</Types>', '<Override PartName="/word/footer1.xml" '
                        'ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.footer+xml"/></Types>')
    files['[Content_Types].xml'] = ct.encode()


def postprocess(path, n_eq):
    zin = zipfile.ZipFile(path)
    files = {i.filename: zin.read(i.filename) for i in zin.infolist()}
    zin.close()
    doc = etree.fromstring(files['word/document.xml'])
    body = doc.find(q('w:body'))
    sect = body.find(q('w:sectPr'))
    # A4, 25 mm margins
    for el in list(sect):
        if el.tag in (q('w:pgSz'), q('w:pgMar')):
            sect.remove(el)
    ps = etree.SubElement(sect, q('w:pgSz')); ps.set(q('w:w'), '11906'); ps.set(q('w:h'), '16838')
    pm = etree.SubElement(sect, q('w:pgMar'))
    for k, v in (('top', '1418'), ('right', '1418'), ('bottom', '1418'), ('left', '1418'),
                 ('header', '709'), ('footer', '709'), ('gutter', '0')):
        pm.set(q('w:' + k), v)
    fr = etree.Element(q('w:footerReference')); fr.set(q('w:type'), 'default'); fr.set(q('{%s}id' % R if False else 'r:id'), 'rIdFooter1')
    sect.insert(0, fr)
    textw = 11906 - 2 * 1418
    bmk = Bookmarks()
    stats = dict(eq=0, ref=0, lnk=0, cap=0)
    # equations
    for p in list(body.iter(q('w:p'))):
        txt = para_text(p).strip()
        m = re.fullmatch(r'\[\[EQ\|(\d+)\|([^\]]*)\]\]', txt)
        if not m:
            continue
        prev = p.getprevious()
        assert prev is not None and prev.find('.//' + q('m:oMathPara')) is not None, 'equation %s not found' % m.group(1)
        parent = prev.getparent()
        idx = parent.index(prev)
        parent.remove(prev)
        parent.insert(idx, eq_table(prev, int(m.group(1)), m.group(2), bmk, textw))
        parent.remove(p)
        stats['eq'] += 1
    # remaining markers
    for p in list(body.iter(q('w:p'))):
        if '[[' in para_text(p):
            split_markers(p, bmk, stats)
    left = [para_text(p) for p in body.iter(q('w:p')) if '[[' in para_text(p)]
    assert not left, left[:3]
    assert stats['eq'] == n_eq, (stats, n_eq)
    # pandoc/texmath schema slips: <m:mcPr> needs count before mcJc; nor and sty are exclusive
    for mc in doc.iter(q('m:mcPr')):
        c = mc.find(q('m:count'))
        if c is not None:
            mc.remove(c); mc.insert(0, c)
    for mr in doc.iter(q('m:rPr')):
        if mr.find(q('m:nor')) is not None:
            for sty in mr.findall(q('m:sty')):
                mr.remove(sty)
    # data tables: content-weighted column widths over the text width; >= 10 columns -> landscape page
    for tbl in list(body.iter(q('w:tbl'))):
        tp = tbl.find(q('w:tblPr'))
        b_ = tp.find(q('w:tblBorders'))
        if b_ is not None and b_.find(q('w:top')).get(q('w:val')) == 'nil':
            continue                                   # equation table
        rows = tbl.findall(q('w:tr'))
        ncol = len(tbl.find(q('w:tblGrid')).findall(q('w:gridCol')))
        land = ncol >= 10
        width = (16838 - 2 * 1418) if land else textw
        wts = []
        for j in range(ncol):
            words, lens = [1], [1]
            for tr in rows:
                tcs = tr.findall(q('w:tc'))
                if j < len(tcs):
                    txt = ''.join(t.text or '' for t in tcs[j].iter(q('w:t'))) + ''.join(
                        t.text or '' for t in tcs[j].iter(q('m:t')))
                    words += [len(w) for w in txt.split()]
                    lens.append(len(txt))
            wts.append(max(max(words), 0.45 * (sum(lens) / len(lens))) + 2)
        cols = [int(width * w / sum(wts)) for w in wts]
        grid = tbl.find(q('w:tblGrid'))
        for g, c in zip(grid.findall(q('w:gridCol')), cols):
            g.set(q('w:w'), str(c))
        for tr in rows:
            for tc, c in zip(tr.findall(q('w:tc')), cols):
                pr = tc.find(q('w:tcPr'))
                if pr is None:
                    pr = etree.Element(q('w:tcPr')); tc.insert(0, pr)
                cw = pr.find(q('w:tcW'))
                if cw is None:
                    cw = etree.Element(q('w:tcW')); pr.insert(0, cw)
                cw.set(q('w:w'), str(c)); cw.set(q('w:type'), 'dxa')
        tw = tp.find(q('w:tblW'))
        if tw is None:
            tw = etree.SubElement(tp, q('w:tblW'))
        tw.set(q('w:w'), str(width)); tw.set(q('w:type'), 'dxa')
        lay = tp.find(q('w:tblLayout'))
        if lay is None:
            lay = etree.Element(q('w:tblLayout'))
            tw.addnext(lay)
        lay.set(q('w:type'), 'fixed')
        if land:
            cap = tbl.getprevious()           # caption paragraph precedes the table
            first = cap if (cap is not None and 'Table' in para_text(cap)[:8]) else tbl
            def sect_par(orient):
                p = etree.Element(q('w:p')); ppr = etree.SubElement(p, q('w:pPr'))
                sp = copy.deepcopy(sect); ppr.append(sp)
                for fr_ in sp.findall(q('w:footerReference')):
                    pass
                pg = sp.find(q('w:pgSz'))
                if orient == 'land':
                    pg.set(q('w:w'), '16838'); pg.set(q('w:h'), '11906'); pg.set(q('w:orient'), 'landscape')
                return p
            first.addprevious(sect_par('port'))     # ends the preceding portrait section
            tbl.addnext(sect_par('land'))           # ends the landscape section holding the table
    files['word/document.xml'] = etree.tostring(doc, xml_declaration=True, encoding='UTF-8', standalone=True)
    add_footer(None, files)
    ct = files['[Content_Types].xml'].decode()
    if 'Extension="png"' not in ct:
        ct = ct.replace('<Default Extension="xml"', '<Default Extension="png" ContentType="image/png" /><Default Extension="xml"')
    files['[Content_Types].xml'] = ct.encode()
    zout = zipfile.ZipFile(path, 'w', zipfile.ZIP_DEFLATED)
    for k, v in files.items():
        zout.writestr(k, v)
    zout.close()
    return stats


def main():
    flat, n_eq = convert_tex()
    ref = reference_doc()
    cmd = ['pandoc', flat, '-f', 'latex', '-t', 'docx', '-o', OUT, '--citeproc',
           '--bibliography', os.path.join(ISA, 'refs.bib'), '--csl', os.path.join(HERE, 'elsevier-harvard.csl'),
           '--reference-doc', ref, '--number-sections', '-M', 'link-citations=true',
           '-M', 'reference-section-title=References', '--resource-path', HERE]
    subprocess.run(cmd, check=True, cwd=HERE)
    stats = postprocess(OUT, n_eq)
    print('written', OUT, stats)


if __name__ == '__main__':
    main()
