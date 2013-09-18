# The following is a tedious way to make a bracket
def bracket(ax, x1, x2, y, em, text):
  xmiddle = 0.5*(x1+x2)
  xleft = 0.5*(x1+xmiddle)
  xight = 0.5*(xmiddle+x2)
  bracket_height = 0.6*em
  ax.text(xmiddle, y-(em+bracket_height), text, horizontalalignment='center',
          transform=ax.transAxes)
  shrink = 0
  rad = 2
  ax.annotate("",
              xy=(x1,y), xycoords='axes fraction',
              xytext=(xleft,y-0.5*bracket_height), textcoords='axes fraction',
              arrowprops=dict(arrowstyle="-",
                               color="0.0",
                               shrinkA=0, shrinkB=shrink,
                               patchA=None,
                               patchB=None,
                               connectionstyle="angle,angleA=0,angleB=-45,rad=%g" % rad,
                               ),
               )
  ax.annotate("",
              xy=(xleft,y-0.5*bracket_height), xycoords='axes fraction',
              xytext=(xmiddle,y-bracket_height), textcoords='axes fraction',
              arrowprops=dict(arrowstyle="-",
                               color="0.0",
                               shrinkA=shrink, shrinkB=0,
                               patchA=None,
                               patchB=None,
                               connectionstyle="angle,angleA=-45,angleB=0,rad=%g" % rad,
                               ),
               )
  ax.annotate("",
              xy=(x2,y), xycoords='axes fraction',
              xytext=(xight,y-0.5*bracket_height), textcoords='axes fraction',
              arrowprops=dict(arrowstyle="-",
                               color="0.0",
                               shrinkA=0, shrinkB=shrink,
                               patchA=None,
                               patchB=None,
                               connectionstyle="angle,angleA=180,angleB=45,rad=%g" % rad,
                               ),
               )
  ax.annotate("",
              xy=(xight,y-0.5*bracket_height), xycoords='axes fraction',
              xytext=(xmiddle,y-bracket_height), textcoords='axes fraction',
              arrowprops=dict(arrowstyle="-",
                               color="0.0",
                               shrinkA=shrink, shrinkB=0,
                               patchA=None,
                               patchB=None,
                               connectionstyle="angle,angleA=45,angleB=0,rad=%g" % rad,
                               ),
               )

