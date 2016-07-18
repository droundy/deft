#!/usr/bin/python

from __future__ import division
from math import pi       # REALLY don't need all of math
import os, numpy as np, sys

#lists all valid scrunched folders in the following format:
#validFolders[m] is the mth base folder (ranked by L)
#validFolders[m][i] is ranked by the order (i=0 most likely means zeroth order)
#validFolders[m][i][n] is ranked by atoms in a box (n=0 does not mean zero atoms)

baseDir=os.getcwd()
def findFolders():
	baseFolders=os.listdir(baseDir)
	global validFoldersss
	validFoldersss=[]
	for folder in baseFolders:
		folder=baseDir+"/"+folder
		if "scrunched-ww" in folder:
			#print folder
			subFolders=os.listdir(folder)
			validFolderss=[]
			for subFolder in subFolders:#expected i=0 i=1 i=... 
				foundL=False
				subFolder=folder+"/"+subFolder
				subNfolders=os.listdir(subFolder)
				validFolders=[]
				for subNfolder in subNfolders:#expected N=3 N=4 N=...
					subNfolder=subFolder+"/"+subNfolder+"/"
					validFolders.append(subNfolder)
				validFolderss.append(validFolders)
			validFoldersss.append(validFolderss)
		sortDone=False
		while sortDone==False:
			sortDone=True
			for i in range(0,len(validFoldersss)-1):
				Lprev=float(validFoldersss[i][0][0].split("-L")[1].split("/i")[0])
				Lnext=float(validFoldersss[i+1][0][0].split("-L")[1].split("/i")[0])
				if Lnext<Lprev:
					sortDone=False
					validFoldersss[i],validFoldersss[i+1]=validFoldersss[i+1],validFoldersss[i]
		
		for i in range(0,len(validFoldersss)):
			sortDone=False
			while sortDone==False:
				sortDone=True
				for j in range(0,len(validFoldersss[i])-1):
					iprev=float(validFoldersss[i][j][0].split("/i")[1].split("/N")[0])
					inext=float(validFoldersss[i][j+1][0].split("/i")[1].split("/N")[0])
					if inext<iprev:
						sortDone=False
						validFoldersss[i][j],validFoldersss[i][j+1]=validFoldersss[i][j+1],validFoldersss[i][j]
		
		for i in range(0,len(validFoldersss)):
			for j in range(0,len(validFoldersss[i])):
				sortDone=False
				while sortDone==False:
					sortDone=True
					for k in range(0,len(validFoldersss[i][j])-1):
						nPrev=float(validFoldersss[i][j][k].split("/N")[1].split("/")[0])
						nNext=float(validFoldersss[i][j][k+1].split("/N")[1].split("/")[0])
						if nNext<nPrev:
							sortDone=False
							validFoldersss[i][j][k],validFoldersss[i][j][k+1]=validFoldersss[i][j][k+1],validFoldersss[i][j][k]
	
