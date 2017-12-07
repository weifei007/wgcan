#!/usr/bin/env python
#-*-coding:utf8-*-

import sys
import re
import pandas as pd

class Hub():
	def __init__(self,model,connect):
		self.model=model
		self.connect=connect
	def _readcon(self,per=5):
		conf=self.connect
		gene_con={}
		with open(conf,'r') as con:
			for line in con:
				if line.startswith("GeneID"):
					continue
				lines=line.strip().split()
				gene_con[lines[0]]=float(lines[1])
		newgene_con=sorted(gene_con.iteritems(),key=lambda a:a[1],reverse=True)
		allnum=len(newgene_con)
		hubs=[]
		n=1
		while n < allnum*per/100:
			nn=n-1
			hubs.append(newgene_con[nn][0])
		return hubs
	def _readmod(self):
		ids=[]
		mf=self.model
		with open(mf,'r') as con:
			for line in con:
				if line.startswith("fromNode"):
					continue
				lines=line.strip().split()
				ids.extend([lines[0],lines[1]])
		idsa=set(ids)
		return idsa
	def _gethub(self,idl,per=5):
		conf=self.connect
		data=pd.read_table(conf,sep="\t",header=0,index_col=0)
		datasel=data.loc[idl,:]
		datasort=datasel.sort_values("kWithin",ascending=False)
		hang=datasort.shape[0]
		getnum=int(hang*per/100)
		id=list(datasort.index[0:getnum])
		ida="\n".join(id)
		return id
	def main(self,weight=0.3):
		modf=self.model
		midl=self._readmod()
		hublist=self._gethub(midl)
		all="Id\tLabel\tcolor\tannot\n"
		idlist=[]
		modcon="Source\tTarget\tWeight\ttype\tfromAltName\ttoAltName\n"
		with open(modf,'r') as con:
			for line in con:
				lines=line.strip().split()
				if line.startswith("fromNode") or float(lines[2])<weight:
					continue
				lb1=""
				c1="0"
				if lines[0] in hublist:
					c1="1"
					lb1=lines[0]
				c2="0"
				lb2=""
				if lines[1] in hublist:
					c2="1"
					lb2=lines[1]
				if lines[0] in hublist: 
					modcon+=line
				else:
					newcon="%s\t%s\t%s\t%s\t%s\t%s\n" % (lines[1],lines[0],lines[2],lines[3],lines[4],lines[5])
					modcon+=newcon	
				if lines[0] in hublist or lines[1] in hublist:
					if lines[0] not in idlist:
						all+=lines[0]+"\t"+lb1+"\t"+c1+"\n"
						idlist.append(lines[0])
					if lines[1] not in idlist:
						all+=lines[1]+"\t"+lb2+"\t"+c2+"\n"
						idlist.append(lines[1])
		na=re.findall(r'cyt-edge_([^\.]+).txt',modf)
		modout=na[0]+"-"+str(weight)+"-edge.csv"
		mout=open(modout,'w')
		mout.write(modcon)
		mout.close()
		out=na[0]+"-"+str(weight)+".id.csv"
		outf=open(out,'w')
		outf.write(all)
		outf.close()
if __name__=="__main__":
	if len(sys.argv) < 3:
		print "\tpython gethubgene.py <cyt-edge_*.txt> <intramodularConnectivity.xls> \n"
		sys.exit()
	All=Hub(sys.argv[1],sys.argv[2])
	All.main(0.3)