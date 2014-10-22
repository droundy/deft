

import markdown

class MathJaxPattern(markdown.inlinepatterns.Pattern):

    def __init__(self):
        markdown.inlinepatterns.Pattern.__init__(self, r'(?<!\\)(\$\$?)(.+?)\2')

    def handleMatch(self, m):
        node = markdown.util.etree.Element('mathjax')
        node.text = markdown.util.AtomicString(m.group(2) + m.group(3) + m.group(2))
        return node

class MathJaxExtension(markdown.Extension):
    def extendMarkdown(self, md, md_globals):
        # Needs to come before escape matching because \ is pretty important in LaTeX
        md.inlinePatterns.add('mathjax', MathJaxPattern(), '<escape')

# for compatibility with python markdown version <= 2.3.1, makeExtension must accept a single argument,
# but for compatibility with python markdown version >= 2.5, MathJaxExtension cannot be passed any arguments.
# as a compromise, makeExtension accepts an optional argument which is never used.
def makeExtension(content=None):
    return MathJaxExtension()


