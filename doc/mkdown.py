
import re, string, os, glob
import markdown as mmdd

import xml.etree.ElementTree as ET

def mkdown(mdfile):
    htfile = mdfile[:-2]+'html'

    f = open(mdfile, 'r')
    mkstr = f.read()
    f.close()

    f = open('sidebar.md', 'r')
    sidebar = f.read()
    f.close()

    f = open('template.html', 'r')
    templatestr = f.read()
    f.close()

    titlere = re.compile(r"^\s*#\s*([^\n]*)(.*)", re.DOTALL)
    title = titlere.findall(mkstr)
    if len(title) == 0:
        title = "Recipes"
    else:
        mkstr = title[0][1]
        title = mmdd.markdown(title[0][0])
        title = title[3:len(title)-4]
    pagetitle = title

    template = string.Template(templatestr)

    f = open(htfile, 'w')
    myhtml = string.replace(template.safe_substitute(
            title = title,
            pagetitle = pagetitle,
            content = mmdd.markdown(mkstr, extensions=['def_list', 'markdown.extensions.toc']),
            sidebar = mmdd.markdown(sidebar, extensions=['def_list'])),
                            '<li><a href="'+mdfile[:-3],
                            '<li><a class="current" href="'+mdfile[:-3])

    ff = open('temp.html', 'w')
    ff.write(myhtml)
    ff.close()
    # etree = ET.fromstring(myhtml)
    # for main in etree.iter('article'):
    #     lastheader = None
    #     for p in list(main):
    #         if p.tag == 'h2' or lastheader is None:
    #             lastheader = ET.SubElement(main, 'div', {'class': 'indivisible'})
    #             lastheader.append(p)
    #             main.remove(p)
    #         elif lastheader is not None:
    #             lastheader.append(p)
    #             main.remove(p)

    # f.write(ET.tostring(etree, encoding='utf-8', method='html'))

    f.write(myhtml)
    f.close()

for f in glob.glob('*.md'):
    mkdown(f)
