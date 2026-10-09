import sys,zipfile,re,html
for p in sys.argv[2:]:
    z=zipfile.ZipFile(p); x=z.read('word/document.xml').decode('utf8')
    out=[]
    for para in re.findall(r'<w:p[ >].*?</w:p>',x,flags=re.S):
        st=re.search(r'<w:pStyle w:val="([^"]+)"',para)
        t=''.join(html.unescape(s) for s in re.findall(r'<w:t[^>]*>([^<]*)</w:t>',para))
        if '<m:oMath' in para: t+=' [MATH:'+''.join(re.findall(r'<m:t[^>]*>([^<]*)</m:t>',para))+']'
        if t.strip(): out.append((('#'+st.group(1)+'# ') if st and 'Head' in st.group(1) else '')+t)
    n=sys.argv[1]+'/'+p.split('/')[-1].replace('.docx','.txt')
    open(n,'w').write('\n'.join(out)); print(n,len(out),sum(len(o) for o in out), 'tables' if '<w:tbl>' in x else '')
