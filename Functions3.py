import os
from Bio import Phylo
from Bio.Phylo.Consensus import *
from io import StringIO
import numpy as np
import numpy
import glob
from scipy import stats
import pandas as pd
from scipy.stats.distributions import chi2
from scipy.stats import fisher_exact
def ReadFasSeq(pos72):
 StLs=[]
 St2Seq={}
 pos72=open(pos72,'r').readlines()
 for i in pos72:
    if i[0]=='>':
       ID=i.strip()
       StLs.append(ID)
       St2Seq[ID]=''
    else: St2Seq[ID]+=i.strip()
 return StLs,St2Seq
def ReadMegSeq(Meg): #input is mega alignment file. out is name2seq dictionary and mega head
  Meg=open(Meg,'r').readlines()
  Read='s'
  out2=''
  NameOrder=[]
  Name2Seq={}
  SeqNum=0
  for i in Meg:
    if i[0]=='#' and i.strip()!='#MEGA' and i.strip()!='#mega' :
      	
        Read='n'
        Name=i.strip()
        NameOrder.append(Name)
        Name2Seq[Name]=''
        SeqNum+=1
    elif Read=='n': Name2Seq[Name]+=i.strip()
    elif Read=='s': out2+=i
  return NameOrder, Name2Seq  
def ReadGV(GV):
    Dec2Anc={}
    NodeMut={}
    Node2In={}
    Edge2In={}	
    GV=open(GV,'r').readlines()#[1:]
    Read='N'
    for i in GV:
        i=i.replace('digraph D {','')	
        if i.find('->')!=-1:
             i0=i.replace(';','').split('->')
             Dec2Anc[i0[1].split('[')[0].strip()]=i0[0].strip()
             Edge2In[i0[0].strip()+'->'+i0[1].split('[')[0].strip()]=i.strip()	
        elif i[0]=='}': Read='D'			
        elif i.strip()!='' and i.find('[')!=-1:
             i0=i.split(' [')
            # print (i0)			 
             Node=i0[0].strip()
             if i0[1].find('\\nMut:')!=-1:			 
                 MutLs=i0[1].split('\\nMut:')[-1].replace('\\n\"]\n','').split(';')
             else: MutLs=[]				 
             MutIn=[]
             for M in MutLs:
                 M=M.strip().split('\\n')			   
                 MutIn+=M	
             NodeMut[Node]=MutIn	
             Node2In[Node]=i		 

	 
        else: pass		
	
    return Dec2Anc,NodeMut,Node2In,Edge2In	
def ReadGV1(GV):
    Dec2Anc={}
    NodeMut={}
    Node2In={}
    Edge2In={}	
    GV=open(GV,'r').readlines()#[1:]
    Read='N'
    for i in GV:
        i=i.replace('digraph D {','')	
        if i.find('->')!=-1:
             i0=i.replace(';','').split('->')
             Dec2Anc[i0[1].split('[')[0].strip()]=i0[0].strip()
             Edge2In[i0[0].strip()+'->'+i0[1].split('[')[0].strip()]=i.strip()	
        elif i[0]=='}': Read='D'			
        elif i.strip()!='' and i.find('[')!=-1:
             i0=i.split(' [')			 
             Node=i0[0].strip()
             MutIn0=i0[1].replace('\"]\n','').split('\\n')[1:]
             MutIn=[]
             for M in MutIn0:
               if M.strip()!='':
                 if M[0]=='B': M='T'+M.replace('Back','').replace(' ','')+'A'
                 elif M[0]=='R': M='Rec A'+M.replace('Rec','').replace(' ','')+'T'
                 else: M='A'+M+'T'
                 MutIn.append(M)
             NodeMut[Node]=MutIn	
             Node2In[Node]=i		 

	 
        else: pass		
	
    return Dec2Anc,NodeMut,Node2In,Edge2In	  
def ReadGVmoa(GV):
    Dec2Anc={}
    NodeMut={}
    Node2In={}
    Edge2In={}	

    GV=open(GV,'r').readlines()#[1:]
    Read='N'
    for i in GV:
        i=i.replace('digraph G {','')	
        if i.find('->')!=-1:
           if i.find('color=\"#ff00005f\"')==-1:
             i0=i.replace(';','').split('->')
             Dec2Anc[i0[1].split('[')[0].strip()]=i0[0].strip()
             Edge2In[i0[0].strip()+'->'+i0[1].split('[')[0].strip()]=i.strip()	
        elif i[0]=='}': Read='D'			
        elif i.strip()!='' and i.find('[')!=-1:
             i0=i.split(' [')			 
             Node=i0[0].strip()
             MutIn=[i0[1].split('\\n')[0].split('label=')[1].replace('\"','').split(',')[0].strip()]
             NodeMut[Node]=MutIn	
             Node2In[Node]=i		 
        else: pass		
	
    return Dec2Anc,NodeMut,Node2In,Edge2In
def ReadNwkWithNodeID(Tree):
    tree = Phylo.read(Tree, "newick")
    C=1
    Name2NodeID={}
    for i in tree.find_clades():
       if i.name==None:
            i.name=C
       Name2NodeID[i.name]=i.comment     
       C+=1
    Tips=tree.get_terminals()
    Dec2Anc={}
    for Tip in Tips:   
       Root2TipLs=tree.get_path(Tip.name)
       D=1
       Len=len(Root2TipLs)
       if Len==1:
           Dec2Anc[Root2TipLs[0].name]='root'       
       while D<Len:
           Dec2Anc[Root2TipLs[D].name]=Root2TipLs[D-1].name
           D+=1
    Dec2Anc[Root2TipLs[0].name]='root'            
    return Dec2Anc,Name2NodeID,Tips	    	
def ReadNodeMap(NodeMap):
    NodeMap=open(NodeMap,'r').readlines()[1:]
    D2A={}
    A2D={}
    N2C={}
    C2N={}
    for i in NodeMap:
        Ls=i.strip().split('\t')
        Line=[]
        for Item in Ls:
            Item=Item.strip().replace(' ','_')
            if Item!='':Line.append(Item)
        N=Line[0]
        C=Line[1]
        N2C[N]=C
        C2N[C]=N
        if Line[3]!='-':
                A2D[C]=[Line[2],Line[3]]
                D2A[Line[2]]=C
                D2A[Line[3]]=C
    return D2A,A2D,N2C,C2N
def ReadCOI(COI):
    COI=open(COI,'r').readlines()
    ParOr=COI[0].strip().split('\t')[1:]
    Len=len(ParOr)
    COI=COI[1:]
    Dec2Anc2COI={}
    for i in COI:
        i=i.split('\t')
        Chi=i[0]
        Dec2Anc2COI[Chi]={}	
        COIs=i[1:]
        c=0
        while c<Len:
           Dec2Anc2COI[Chi][ParOr[c]]=float(COIs[c])
           c+=1
    return Dec2Anc2COI	
def ReadFreTa(Ta):
    Ta=open(Ta,'r').readlines()[1:]
    Pos2Fre={}
    for i in Ta:
        i=i.split('\t')
        Pos2Fre[int(i[0])]=float(i[4])
    return Pos2Fre
def ListCol(File):
  File=open(File,'r').readlines()
  NameOrder,Name2Col=GetHead(File[0])
  File=File[1:]
  Tu2Freq={}
  for Tu in NameOrder:
    Tu2Freq[Tu]=[]
  for i in File:
    i=i.strip().split('\t')
    for Tu in Name2Col:
        Tu2Freq[Tu].append(i[Name2Col[Tu]])
  return Tu2Freq
def GetHead(Head):
    Head=Head.strip().split('\t')
    Len=len(Head)
    c=0
    Name2Col={}
    NameOrder=[]	
    while c<Len:
        Name2Col[Head[c]]=c
        NameOrder.append(Head[c])		
        c+=1
    return NameOrder,Name2Col	
def MOA2GV(MOA):
    GV=MOA[:-4]+'.gv'
    Dec2Anc,NodeInf,Node2In,Edge2In=ReadGV(MOA)
    ID2Lab={}
    for N in Node2In:
       Lab=Node2In[N].split('label=')[-1].split(',')[0].split('\\n')[0]
       ID2Lab[N]=Lab.replace('\"','')	
    out='digraph G {\n'
    for N in Node2In:
       In=Node2In[N]	
       out+=ID2Lab[N]+In[len(N):]
    for E in Edge2In:
      if Edge2In[E].find('label=')!=-1 or E.find('root')!=-1:	
       E=E.split('->')
       out+=ID2Lab[E[0].strip()]+'->'+ID2Lab[E[1].replace(';','').strip()]+';\n'
    out+='}'
    GetOut(GV,out)	
def GetVAFpos(seq_annoFile,VAF):	#Hap+os.sep+'sequence_annotations_Force.txt'
    Pos2Cou={}
    File=open(seq_annoFile,'r').readlines()[1:]
    for i in File:
        i=i.split('\t')
        P=i[0]
        if P not in Pos2Cou: Pos2Cou[P]={'Tot':0,'Var':{}}
        Ref=i[1]
        Alt=i[2]
        Cou=int(i[3])		
        if Ref!=Alt and Ref!='?' and Alt!='?' and Ref!='-' and Alt!='-':
             Pos2Cou[P]['Var'][Alt]=Cou 
        Pos2Cou[P]['Tot']+=Cou
    Pls=[]		
    for P in Pos2Cou:
        Var2Cou=Pos2Cou[P]['Var']
        Tot=Pos2Cou[P]['Tot']		
        for	V in Var2Cou:
            Pro=1.0*Var2Cou[V]/Tot	
            if Pro>VAF: 
                Pls.append(int(P))			
    Pls=list(set(Pls))
    Pls.sort()	
    print (len(Pls))	
    GetOut(seq_annoFile[:-4]+'_Force.txt','\n'.join(map(str,(Pls))))	
		
def GetMostSim(Hap,RefHap):
    DifC2HapLs={}
    Len=len(Hap)
    for Ref in RefHap:
        RefSeq=RefHap[Ref]
        if len(RefSeq)!=Len:
            print (Hap,RefSeq,Len,len(RefSeq))
            open('a','r').readlines()
        else:
           DC=CountDifNum(Hap,RefSeq)
           DifC2HapLs[DC]=DifC2HapLs.get(DC,[])+[Ref.replace('>','')]
    DifCLs=list(DifC2HapLs.keys())
    DifCLs.sort()
    DifC=DifCLs[0]
    Gid=';'.join(DifC2HapLs[DifC])

    return Gid,DifC
def CountDifNum_excMiss(Seq0,Seq1):
            Len=len(Seq0)		
            Dif=0
            c=0
            while c<Len:
                if Seq0[c]!='?' and Seq1[c]!='?' and Seq0[c]!='-' and Seq1[c]!='-':			
                    if Seq0[c]!=Seq1[c]: Dif+=1
                c+=1
            return Dif		
def CountDifNum(Seq0,Seq1):
            Len=len(Seq0)
            Dif=0
            c=0
            while c<Len:
                if Seq0[c]!=Seq1[c]: Dif+=1
                c+=1
            return Dif
def CountDifPos_from1(Seq0,Seq1):
            Len=len(Seq0)
            Dif=[]
            c=0
            while c<Len:
                if Seq0[c]!=Seq1[c]: Dif.append(c+1)
                c+=1
            return Dif
def makeRaxMLin(OriDic,Qseq,OutFas,Add):
    for Q in Qseq:
        InFas='>Query_'+Q.replace('>','')+'\n'+Qseq[Q]+Add+'\n'
    for O in OriDic:
      InFas+=O.replace('__','_')+'\n'+OriDic[O]+Add+'\n'
    GetOut(OutFas,InFas)

def cleanRaxMLtree(Tree,Out,OutG):
    Len=len(Tree)
    c=0
    Read='s'
    New=''
    while c<Len:
       if Tree[c]=='[':
           Read='r'
       elif Read=='r':
          if Tree[c]==']':
           Read='s'
       elif Read=='s': New+=Tree[c]
       c+=1
    GetOut('Unroot.nwk',New)
    root_tree('Unroot.nwk',OutG)
    GetOut(Out,open('Unroot_rooted.nwk','r').readlines()[0])
    os.remove('Unroot.nwk')
def prune_tree(OriNwk,ExtraLs,OutSeq):
   if os.path.exists('rooted.nwk')==True: os.remove('rooted.nwk')
   trees = list(Phylo.parse(OriNwk, 'newick'))
   for tree in trees:
       tree = tree.root_with_outgroup({'name': OutSeq})
   Phylo.write(trees, 'rooted.nwk', "newick")
   TreeLs=open('rooted.nwk','r').readlines()
   Cou=1
   for Tst in TreeLs:
    tree=Phylo.read(StringIO(Tst.strip()), "newick")
    for tip in ExtraLs:
         tree.prune(tip)
    Phylo.write(tree, OriNwk[:-4]+'_'+str(Cou)+'_prune.nwk','newick')
    Cou+=1

def root_tree(OriNwk,Root):
    Out=OriNwk[:-4]+'_rooted.nwk'
    trees = list(Phylo.parse(OriNwk, 'newick'))
    for tree in trees:
       tree = tree.root_with_outgroup({'name': Root})
    Phylo.write(trees, Out, "newick")
def Getalldec(N,Dec2Anc):
    Anc2Dec=GetAnc2Dec(Dec2Anc)
    DecLs=Anc2Dec.get(N,[])
    DecLsAll=Anc2Dec.get(N,[])	
    Added='y'
    while DecLs!=[]:
        NewDecLs=[]
        for D in DecLs:
            NewDecLs+=Anc2Dec.get(D,[])
 
        DecLsAll+=NewDecLs
        DecLs=NewDecLs
    return DecLsAll	
def GetAnc2Dec(Dec2Anc):
    Anc2Dec={}
    for D in Dec2Anc:
        A=Dec2Anc[D]	 
        Anc2Dec[A]=Anc2Dec.get(A,[])+[D]	
    return Anc2Dec	
def GetTipLs(Dec2Anc):
    Anc2Dec=InvertDic1(Dec2Anc)
    TipLs=[]
    for Dec in Dec2Anc:
        if Dec not in Anc2Dec: TipLs.append(Dec)
    return TipLs        
def ClassifyDif(AncSeq,DecSeq,PosIDOrder):
    ForLs=[]
    BackLs=[]
    Len=len(AncSeq)
    c=0
    while c<Len:
        if AncSeq[c]=='A' and DecSeq[c]=='T': ForLs.append(PosIDOrder[c].strip())
        elif AncSeq[c]=='T' and DecSeq[c]=='A': BackLs.append(PosIDOrder[c].strip())
        c+=1
    return ForLs,BackLs 
def ClassifyDif1(AncSeq,DecSeq,PosIDOrder):
    ForLs=[]
    BackLs=[]
    Len=len(AncSeq)
    print (Len,len(DecSeq))
    c=0
    while c<Len:
        if AncSeq[c] != DecSeq[c]: ForLs.append(AncSeq[c]+PosIDOrder[c].strip()+DecSeq[c])	
        c+=1
    return ForLs,BackLs 	
def MutTree2dot(File):
   Out=File[:-4]+'.gv'

   File=open(File,'r').readlines()[1:]
   Nodes=''
   Paths=''
   for i in File:
       i=i.split('\t')
       
       if len(i)>5: Nodes+=i[0]+	' [label=\"'+i[0]+'\\nMut:'+i[2]+'\\nBack:'+i[3]+'\\nRec:'+i[4]+'\\nCellC:'+i[5].strip()+'\"]\n'	
       else: Nodes+=i[0]+	' [label=\"'+i[0]+'\\nMut:'+i[2]+'\\nBack:'+i[3]+'\\nRec:'+i[4].strip()+'\"]\n'	   
       Paths+=i[1]+'->'+i[0]+'\n'	   

   out='digraph D {\n'+Nodes+Paths+'}\n'
   GetOut(Out,out)
def MutTree2dot1(File):
   Out=File[:-4]+'.gv'

   File=open(File,'r').readlines()[1:]
   Nodes=''
   Paths=''
   for i in File:
       i=i.replace('.','').split('\t')
       
       #if len(i)>5: Nodes+=i[0]+	' [label=\"'+i[0]+'\\nMut:'+i[2]+'\\nBack:'+i[3]+'\\nRec:'+i[4]+'\\nCellC:'+i[5].strip()+'\"]\n'	
       Nodes+=i[0]+	' [label=\"'+i[0]+'\\nMut:'
       MutLs=i[2].split(';')
       In=[]
       c=0
       for Mut in MutLs:
           if c>10: 
               Nodes+=';'.join(In)+'\\n'
               In=[]
               c=0
           
           In.append(Mut)
           c+=1		   
       if In!=[]:    Nodes+=';'.join(In)+'\\n'               		   
       Nodes+='\"]\n'	   
       Paths+=i[1]+'->'+i[0]+'\n'	   

   out='digraph D {\n'+Nodes+Paths+'}\n'
   GetOut(Out,out)
def GetComAnc(NodeLs,Dec2Anc):
    Node2C={}
    for N in NodeLs:
        Node2C[N]=Node2C.get(N,0)+1
        while N in Dec2Anc:
            N=Dec2Anc[N]
            Node2C[N]=Node2C.get(N,0)+1
    Com=[]
    for N in Node2C:			
        if Node2C[N]>1: Com.append(N)
    return Com
def getCOI(TarM,AncMLs,Dec2Anc2COI):
    AncCOI=[]
    for A in AncMLs:
      if A!='':	
       if TarM.replace('Rec','')[1:]!='' and A.replace('Rec','')[1:]!='':	
        if TarM.replace('Rec','')[1:] in Dec2Anc2COI:
           COI0=Dec2Anc2COI[TarM.replace('Rec','')[1:]]
           if A.replace('Rec','')[1:] in COI0:           
               COI=Dec2Anc2COI[TarM.replace('Rec','')[1:]][A.replace('Rec','')[1:]]
               AncCOI.append(A+'='+str(COI))
    return AncCOI	
def PruneTip(GV,Fas,Nwk,Anno):
    Dec2Anc,NodeMut,Node2In,Edge2In	=ReadGV(GV)
    KeepTH=[]
    PruneTH=[]
    PruneTHNwk=[]
    out='digraph D {\n'
    for TH in NodeMut:
      if TH!='root':
        MutLs=NodeMut[TH]
        if MutLs!=[''] and len(MutLs)>0: 
            KeepTH.append(TH)
            out+=Node2In[TH]
        else: 
           if TH[:5]!='Node_': PruneTHNwk.append(TH)
           PruneTH.append(TH)       
      else: out+=Node2In[TH]
    for E in Edge2In:
      if PruneTH.count(E.split('->')[1].strip())==0:
        out+=Edge2In[E]+'\n'
    out+='}\n'    
    GetOut(GV[:-3]+'1.gv',''.join(open(GV,'r').readlines()))    
    KeepTH.append('Outgroup')
    print (PruneTH)
    print ('prune fasta')
    THls,TH2Seq=ReadFasSeq(Fas)
    out=''
    for TH in THls:
        if PruneTH.count(TH.replace('>',''))==0:
            out+=TH+'\n'+TH2Seq[TH]+'\n'
        else: print ('prune',TH)
    GetOut(Fas[:-6]+'COI.fasta',out)
    print ('prune tree')
    prune_tree(Nwk,PruneTHNwk,'Outgroup')
    print ('prune annotation')    
    Anno=open(Anno,'r').readlines()
    out=Anno[0]
    Anno=Anno[1:]
    for i in Anno:
        TH=i.split('\t')[2].strip()
        if PruneTH.count(TH)==0: out+=i        
    GetOut(Nwk.replace(Nwk.split(os.sep)[-1],'')+'TopHapCOI_anno.txt',out)#'TopHapNA__2022_1_FillAllCell_anno1_COI.txt')    
def CleanMutTree(GV,MutIDLs,Nwk):
    OutDir=GV.replace(GV.split(os.sep)[-1],'')
    Out=OutDir+'MutTree.gv'
    print (Out)
    COI=OutDir+'MOAout'+os.sep+'COI_matrix.txt'
    COIout=COI[:-11]+'.txt'
    CloTa=OutDir+'TopHapCOI_anno.txt'#'TopHapNA__2022_1_FillAllCell_anno1_COI.txt'
    Dec2Anc,NodeMut,Node2In,Edge2In	=ReadGV(GV)
    
    Dec2Anc2COI=ReadCOI(COI)
    
    if os.path.exists(MutIDLs)==True:
        MutIDOrder=open(MutIDLs,'r').readlines()
    else: 
        MutLs0=Dec2Anc2COI.keys()
        print (MutLs0)
        #open('a','r').readlines()
        MutLs=[]
        for i in MutLs0:
    
            MutLs.append(int(i[:-1]))
        MutLs.sort()
        MutIDOrder= list(range(0,MutLs[-1]+1))
        print (MutIDOrder)
    COI=open(COI,'r').readlines()
    Head=COI[0].strip().split('\t')
    out=Head[0]
    Head=Head[1:]
    for i in Head:
        Pos=int(i[:-1])
        out+='\t'+ str(MutIDOrder[Pos]).strip()+'_'+i 
    out+='\n'
    COI=COI[1:]
    for i in COI:
        Pos=int(i.split('\t')[0][:-1])
        out+=str(MutIDOrder[Pos]).strip()+'_'+i 
    GetOut(COIout,out)    
    KeepTipLs=[]
    CloTa=open(CloTa,'r').readlines()[1:]
    for i in CloTa:
        i=i.split('\t')[2]
        KeepTipLs.append(i.strip())
    KeepTipLs=list(set(KeepTipLs))
    print ('clone list',KeepTipLs)    
    Anc2Dec=InvertDic1(Dec2Anc)
    TipLs=GetTipLs(Dec2Anc)
    
    PruneLs=[]
    Mut2C={}
    for N in NodeMut:
        MutLs=NodeMut[N]
        for M in MutLs:
            Mut2C[M]=Mut2C.get(M,0)+1
    for Dec in Dec2Anc:
        if Dec[:3]=='All' and KeepTipLs.count(Dec)==0: PruneLs.append(Dec)
        
    for Dec in Dec2Anc:
        if Dec[:5]=='Node_':
            Lins=Anc2Dec[Dec]
            GoodLn=0
            for Lin in Lins:
                AllDecLs=Getalldec(Lin,Dec2Anc)
                AllDecLs.append(Lin)
                Good='n'
                for D in AllDecLs:
                     if D[:3]=='All' and PruneLs.count(D)==0: Good='y'
                if Good=='y': GoodLn+=1   
            if GoodLn<2: PruneLs.append(Dec)
    print ('prune node',PruneLs)   
    NewNode2In={}
    NewEdge2In={}
    PruneLs=[]
    for Node in Node2In:
       In=Node2In[Node]
       if In.find('Mut:')==-1:NewNode2In[Node]=In
       else:    
        MutLs=In.split('Mut:')[-1].replace('\"]','').split(';')
        NewMutLs=[]
        for M0 in MutLs:
          M0=M0.strip().split('\\n')
          for M in M0:
           if M.strip()!='' and M.strip().replace('Rec','')!='':
            Pos=M.strip().replace('Rec','').strip()[1:][:-1]
            print ( Pos,M)
            MID=str(MutIDOrder[int(Pos)]).strip()
            GoodBack='y'
            if Mut2C[M]>1: MID='Rec '+MID #M.find('Rec')!=-1: MID='Rec '+MID
            if M.strip()[-1]=='A':
                 ForFind='n'
                 Anc=Dec2Anc[Node]
                 AncMutLs=NodeMut[Anc]
                 if AncMutLs.count('A'+Pos+'T')!=0: ForFind='y'
                 while Anc in Dec2Anc:
                     Anc=Dec2Anc[Anc]
                     AncMutLs=NodeMut[Anc]
                     if AncMutLs.count('A'+Pos+'T')!=0: ForFind='y'                       
                 MID='Back '+MID
                 if ForFind=='n': GoodBack='n'
            if GoodBack=='y': NewMutLs.append(MID)#.replace(Pos,':'+MID+':').strip())
            else: 
                print ('remove due to no forward mutation',MID)
        if (NewMutLs==[] or NewMutLs==['']) and Node[:3]=='All': PruneLs.append(Node)         
        NewIn=In.split('Mut:')[0]+'\\n'.join(NewMutLs)+'\"]\n'    
        NewNode2In[Node]=NewIn
    for E in Edge2In:
        In=Edge2In[E]
        COIls=In.split('label=\"')[-1].replace('\"];','').split('\\n')
        NewEdge2In[E]=E+';\n'#NewIn    
    out='digraph D {\n'
    for N in NewNode2In:
        out+=NewNode2In[N]
     
    for E in NewEdge2In:
        out+=NewEdge2In[E]+'\n'
    out+='}'
    GetOut(Out,out)	
    os.system('dot -Tpng '+Out+' -o '+Out[:-3]+'.png')   
    print ('pruen list',PruneLs)    
    prune_tree(Nwk,PruneLs,'Outgroup')     
def CleanMutTree1(GV,MutIDLs,CloTa,SNVc): #CleanMutTree1(TopHapMOAgv,MutIDFile,TopHapMOACellAnno) #output '_final.gv' 
    OutDir=GV.replace(GV.split(os.sep)[-1],'')
    Out=GV[:-3]+'_final.gv'#OutDir+'MutTree.gv'
    print (Out)
    COI=OutDir+'MOAout'+os.sep+'COI_matrix.txt'
    COIout=COI[:-11]+'.txt'
    Dec2Anc,NodeMut,Node2In,Edge2In	=ReadGV1(GV)  
    Dec2Anc2COI=ReadCOI(COI)   
    if os.path.exists(MutIDLs)==True:
        MutIDOrder=open(MutIDLs,'r').readlines()
    else: 
        MutIDOrder= list(range(0,SNVc))
        print (MutIDOrder)
    COI=open(COI,'r').readlines()
    Head=COI[0].strip().split('\t')
    out=Head[0]
    Head=Head[1:]
    for i in Head:
        Pos=int(i[:-1])
        out+='\t'+ str(MutIDOrder[Pos]).strip()+'_'+i 
    out+='\n'
    COI=COI[1:]
    for i in COI:
        Pos=int(i.split('\t')[0][:-1])
        out+=str(MutIDOrder[Pos]).strip()+'_'+i 
    GetOut(COIout,out)    
    KeepTipLs=[]
    CloTa=open(CloTa,'r').readlines()[1:]
    for i in CloTa:
        i=i.split('\t')[1]
        KeepTipLs.append(i.strip())
    KeepTipLs=list(set(KeepTipLs))
    print ('clone list',KeepTipLs)    
    Anc2Dec=InvertDic1(Dec2Anc)
    TipLs=GetTipLs(Dec2Anc)
    
    PruneLs=[]
    Mut2C={}
    for N in NodeMut:
        MutLs=NodeMut[N]
        for M in MutLs:
            Mut2C[M]=Mut2C.get(M,0)+1
    print ('mutID to mut count',Mut2C)        
    for Dec in Dec2Anc:
        if KeepTipLs.count(Dec)==0: PruneLs.append(Dec)
    print ('prune node',PruneLs)   
    NewEdgeIn=[]
    NewNode2Mut={}
    for TN in TipLs:
       if TN not in PruneLs:
           TN1=TN
           NodeID=TN            
           Muts=RenameMut(NodeMut[TN1],Mut2C,MutIDOrder)
           
           while TN1 in Dec2Anc:
              TN1=Dec2Anc[TN1]
              DecC=len(Anc2Dec[TN1])
              if (TN1 in PruneLs or TN1[:3]!='All') and DecC==1:
                  Muts+=RenameMut(NodeMut[TN1],Mut2C,MutIDOrder)
              else: 
                  NewNode2Mut[NodeID]=Muts
                  NewEdgeIn.append(TN1+'->'+NodeID)
                  NodeID=TN1 
                  Muts=RenameMut(NodeMut[TN1],Mut2C,MutIDOrder)
           NewNode2Mut[NodeID]=Muts
    
    out=['digraph D {\n']
    for N in NewNode2Mut:
        if N[:3]!='All': Lab='Node_'+N
        else: Lab=N
        out+=[N+' [label=\"'+Lab+'\\n'+'\\n'.join(NewNode2Mut[N])+'\"]\n']
    NewEdgeIn=list(set(NewEdgeIn)) 
    for E in NewEdgeIn:
        out+=[E+';\n']
    out+=['}']
    GetOut(Out,''.join(out))	 
def RenameMut(MutLs,Mut2C,IDOr):
    NewLs=[]
    for M in MutLs:
        if M[0]=='B' or M[-1]=='A':Back='Back '
        else: Back=''
        if M.find('root')==-1:
           #print (M)
           Mid=int(M.replace('Rec','').replace('Back','').replace('A','').replace('T',''))
           if Mut2C[M]>1: NewLs.append('Rec '+Back+str(IDOr[Mid]).strip())
           else: NewLs.append(Back+str(IDOr[Mid]).strip())
    return NewLs    
def GV2Fas(PosOr0):
    MainFol=PosOr0.replace(PosOr0.split(os.sep)[-1],'')
    COImatrix=MainFol+'TopHapOut'+os.sep+'MOAout'+os.sep+'COI_matrix.txt'
    OutGseq=MainFol+'OutG.fasta' #haplotype
    Nwk=MainFol+'TopHapOut'+os.sep+'TopHap_allMP_AncMP.nwk'
    GV=MainFol+'TopHapOut'+os.sep+'TopHap_MutTree_prune_COI3_ave3_rec_back_COI3_ave3.gv'
    OutF=MainFol+'TopHapOut'+os.sep+'TopHapCOI.fasta'
    CellAnnoTar=MainFol+'TopHapOut'+os.sep+'TopHap*_FillAllCell_anno1.txt'
    PosOr=open(PosOr0,'r').readlines()
    Len=len(PosOr)
    OutGLs,OutG2Seq=ReadFasSeq(OutGseq)
    OutSeq=OutG2Seq[OutGLs[0]]
    print (OutSeq)
    if len(OutSeq)!=Len: open('a','r').readlines()
    
    Dec2Anc,NodeMut,Node2In,Edge2In	=ReadGV(GV)
    print (GV, NodeMut)
    RmTH=[]
    BackRecTH=[]
    out=''
    outRB=''
    for N in Dec2Anc:
        if N.find('Node_')==-1:
            TH='>'+N	
            MLs=NodeMut[N]
            Mut=0
            ForMut=0		
            for M in MLs:
               if M.strip()!='': Mut+=1
               if M.strip()!='': 
                   if M[0]=='T' or M.find('Rec')!=-1: ForMut+=1		   
            if Mut==0: RmTH.append(N)
    		
            else:
               AllMut=MLs
               while N in Dec2Anc:
                   N=Dec2Anc[N]
                   Dec2Anc,NodeMut,Node2In,Edge2In	=ReadGV(GV)
                   if N in NodeMut: 
                       AllMut+=NodeMut[N]
                       print (N,NodeMut[N])
               AllMut.reverse()
               NewSeq=list(OutSeq)
               NewSeqRB=list(OutSeq)	
               print (TH,AllMut)               
               for M in AllMut:
                 if M.strip()!='':		   
                   Pos=M.replace('Rec','')[1:][:-1]
                    			   
                   if PosOr.count(Pos)!=0: P=PosOr.index(Pos)
                   else: P=PosOr.index(Pos+'\n')
                   NewSeq[P]=M[-1]
                   if M[0]=='A': NewSeqRB[P]=M[-1]			   
               NewSeq=''.join(NewSeq)	
               NewSeqRB=''.join(NewSeqRB)		   
               print (NewSeq)
               out+=TH+'\n'+NewSeq+'\n'	   
               if ForMut==Mut: 
    		   
                    BackRecTH.append(TH.replace('>',''))	
               else: outRB+=TH+'\n'+NewSeqRB+'\n'	   
    out+=OutGLs[0]+'\n'+OutSeq+'\n'	 
    outRB +=OutGLs[0]+'\n'+OutSeq+'\n'	 
    GetOut(OutF,out)	
    GetOut(OutF[:-6]+'RB.fasta',outRB)
    print ('prune TopHap haplotype',RmTH)
    if os.path.exists(Nwk)==True:
        prune_tree(Nwk,RmTH,OutGLs[0].replace('>',''))
        if os.path.exists(Nwk[:-17]+'COI.nwk')!=True: os.system('megacc -a analyze_user_tree_MP__nucleotide.mao -d '+OutF+' -t '+Nwk[:-4]+'_1_prune.nwk -o '+Nwk[:-17]+'COI.nwk')
        TopHap_DrawMutTree2(OutF,PosOr0,Nwk[:-17]+'COI.nwk','Outgroup')        
        MapCOI3(COImatrix,OutF[:-6]+'_MutTree.gv')        	
        AveCOI3(OutF[:-6]+'_MutTree_COI3.gv')
    CellAnnoLs=glob.glob(CellAnnoTar)           		   
    print (CellAnnoLs)
    for CellAnno in CellAnnoLs:
        Out=CellAnno[:-4]+'_COI.txt'
        CellAnno=open(CellAnno,'r').readlines()[1:]
        out='Group\tSeq\tTopHapID\n'
        for i in CellAnno:
            if RmTH.count(i.split('\t')[2].strip())==0: out+=i	
    
        GetOut(Out,out)	
    RmTH+=BackRecTH
    RmTH=list(set(RmTH))
def GV2Fas1(GV,PosOr0):
    OutF=GV[:-3]+'.fasta'
    OutP=GV[:-3]+'_Pos.txt'
    if os.path.exists(PosOr0)==True:
        PosOr=open(PosOr0,'r').readlines()
    else: 
        PosOr=range(0,99999)
        PosOr=[str(x) for x in PosOr]
    Len=len(PosOr)
    OutSeq='A'*Len    
    Dec2Anc,NodeMut,Node2In,Edge2In	=ReadGV1(GV)
    print (GV, NodeMut)
    RmTH=[]
    BackRecTH=[]
    out=''
    outRB=''
    TH2Seq={}
    for N in Dec2Anc:
        
        if N.find('Node_')==-1:
            TH='>'+N	
            MLs=NodeMut[N]
            Mut=0
            ForMut=0		
            for M in MLs:
               if M.strip()!='': Mut+=1
               if M.strip()!='': 
                   if M.find('Back')!=-1 or M.find('Rec')!=-1: ForMut+=1		   
            if Mut==0: RmTH.append(N)
    		
            else:
               AllMut=MLs
               while N in Dec2Anc:
                   N=Dec2Anc[N]
                   Dec2Anc,NodeMut,Node2In,Edge2In	=ReadGV1(GV)
                   if N in NodeMut: 
                       AllMut+=NodeMut[N]
                       print (N,NodeMut[N])	
               AllMut.reverse()
               NewSeq=list(OutSeq)
               NewSeqRB=list(OutSeq)	              
               for M in AllMut:
                 if M.strip()!='':		   
                   Pos=M.replace('Rec','').replace('Back','').strip().replace('A','').replace('T','').replace('?','')
                    			   
                   if PosOr.count(Pos)!=0: P=PosOr.index(Pos)
                   else: P=PosOr.index(Pos+'\n')
                   NewSeq[P]='T'
                   if M.find('Back')!=-1: NewSeq[P]='A'			   
               NewSeq=''.join(NewSeq)			   
               out+=TH+'\n'+NewSeq+'\n'
               TH2Seq[TH]=NewSeq	   
    print ('extract variable sites')
    MutP ,ExtSeqDic =	ExtractVarPos(TH2Seq) 
    out='>Outgroup\n'+('A'*len(MutP))+'\n'	
    for TH in ExtSeqDic:
        out+=TH+'\n'+ExtSeqDic[TH]+'\n'
    GetOut(OutF,out)
    out=''
    for M in MutP:
        out+= PosOr[M]+'\n'
    GetOut(OutP,out)
  
    return MutP,ExtSeqDic
def GV2Fas2(GV,PosOr0):
    OutF=GV[:-3]+'.fasta'
    OutP=GV[:-3]+'_Pos.txt'
    if os.path.exists(PosOr0)==True:
        PosOr=open(PosOr0,'r').readlines()
    else: 
        PosOr=range(0,99999)
        PosOr=[str(x) for x in PosOr]
    Len=len(PosOr)
    OutSeq='A'*Len    
    Dec2Anc,NodeMut,Node2In,Edge2In	=ReadGV1(GV)
    print (GV, NodeMut)
    RmTH=[]
    BackRecTH=[]
    out=''
    outRB=''
    TH2Seq={}
    Seq2THls={}
    for N in Dec2Anc:
        
        if N.find('Node_')==-1:
            TH='>'+N	
            MLs=NodeMut[N]
            Mut=0
            ForMut=0		
            for M in MLs:
               if M.strip()!='': Mut+=1
               if M.strip()!='': 
                   if M.find('Back')!=-1 or M.find('Rec')!=-1: ForMut+=1		   
            Go='y' 
            if Go=='y':            
               AllMut=MLs
               while N in Dec2Anc:
                   N=Dec2Anc[N]
                   Dec2Anc,NodeMut,Node2In,Edge2In	=ReadGV1(GV)
                   if N in NodeMut: 
                       AllMut+=NodeMut[N]
                       print (N,NodeMut[N])
               AllMut.reverse()
               NewSeq=list(OutSeq)
               NewSeqRB=list(OutSeq)	             
               for M in AllMut:
                 if M.strip()!='':		   
                   Pos=M.replace('Rec','').replace('Back','').strip().replace('A','').replace('T','').replace('?','')
                    			   
                   if PosOr.count(Pos)!=0: P=PosOr.index(Pos)
                   else: P=PosOr.index(Pos+'\n')
                   NewSeq[P]='T'
                   if M.find('Back')!=-1: NewSeq[P]='A'			   
               NewSeq=''.join(NewSeq)			   
               out+=TH+'\n'+NewSeq+'\n'
               TH2Seq[TH]=NewSeq
               Seq2THls[NewSeq]=Seq2THls.get(NewSeq,[])+[TH]               
    print ('extract variable sites')
    MutP ,ExtSeqDic =	ExtractVarPos(TH2Seq) 
    ExpSeq2THls=InvertDic1(ExtSeqDic)
    out='>Outgroup\n'+('A'*len(MutP))+'\n'	
    for Seq in ExpSeq2THls:
        TH=''.join(ExpSeq2THls[Seq])
        out+='>'+TH.replace('>','')+'\n'+Seq+'\n'
    GetOut(OutF,out)
    out=''
    for M in MutP:
        out+= PosOr[M].strip()+'\n'
    GetOut(OutP,out)
  
    return MutP,ExtSeqDic    
def GV2FasMoa(GV,SNVc):
    OutF=GV[:-3]+'.fasta'
    OutP=GV[:-3]+'_Pos.txt'

    PosOr=range(0,SNVc)
    PosOr=[str(x) for x in PosOr]
    Len=len(PosOr)
    OutSeq='A'*Len    
    Dec2Anc,NodeMut,Node2In,Edge2In	=ReadGV1(GV)
    RmTH=[]
    BackRecTH=[]
    out=''
    outRB=''
    TH2Seq={}
    for N in Dec2Anc:
               TH='>'+N	
               N1=Dec2Anc[N]
               AllMut=[int(N.replace('\"',''))]
               if N1 !='\"root\"' and N1 !='root':AllMut.append(int(N1.replace('\"','')))
               while N1 in Dec2Anc:
                   N1=Dec2Anc[N1]
                   if N1 !='\"root\"' and N1 !='root':AllMut.append(int(N1.replace('\"','')))
               AllMut.reverse()
               NewSeq=list(OutSeq)
               NewSeqRB=list(OutSeq)	             
               for P in AllMut:
                   NewSeq[P]='T'
               NewSeq=''.join(NewSeq)			   
               out+=TH+'\n'+NewSeq+'\n'
               TH2Seq[TH]=NewSeq	   
    print ('extract variable sites')
    MutP ,ExtSeqDic =	ExtractVarPos(TH2Seq) 
    out='>Outgroup\n'+('A'*len(MutP))+'\n'	
    for TH in ExtSeqDic:
        out+=TH+'\n'+ExtSeqDic[TH]+'\n'
    GetOut(OutF,out)

    GetOut(OutP,'\n'.join(map(str,MutP)))
    print (len(MutP))

    return MutP,ExtSeqDic    
def GetMutPos(Seq):
   Len=len(Seq)
   c=0
   Ls=[]
   while c<Len:
      if Seq[c]=='T': Ls.append(c)
      c+=1 
   return Ls       
def ExtractVarPos(CloSeq):
    MutP=[]
    CloneC=0
    for Clo in CloSeq:
        MutP+=GetMutPos(CloSeq[Clo])
        if CloSeq[Clo].find('T')!=-1: CloneC+=1
    MutP=list(set(MutP))
    MutP.sort()
    print (MutP)
    NewSeqDic={}
    for Seq in CloSeq:
        NewSeq=''
        for M in MutP:
            NewSeq+=CloSeq[Seq][M]
        NewSeqDic[Seq]=NewSeq    
    return MutP ,NewSeqDic  
def InvertDic(St2Seq):
 Hap2ID={}
 for St in St2Seq:
    Seq=St2Seq[St]
    Hap2ID[Seq]=St
 return Hap2ID
def InvertDic1(St2Seq):
 Hap2ID={}
 for St in St2Seq:
    Seq=St2Seq[St]
    Hap2ID[Seq]=Hap2ID.get(Seq,[])+[St]
 return Hap2ID
def InvertDic2Ls(St2Seq):
 Hap2ID=[]
 for St in St2Seq:
    Seq=St2Seq[St]
    if Hap2ID.count(Seq)==0:	
        Hap2ID.append(Seq)
 	
 return Hap2ID
def InvertDic3(St2Seq):
 Hap2ID={}
 for St in St2Seq:
    SeqLs=St2Seq[St]
    for Seq in SeqLs:	
        Hap2ID[Seq]=Hap2ID.get(Seq,[])+[St]
 return Hap2ID 
def CountDifNum_RmMiss_AT(Seq0,Seq1): #exp, obs
            Len=len(Seq0)
         		
            Dif={'TP':0,'FN':0,'FP':0,'TN':0,'Tot':0}
            c=0
            if len(Seq0)!=len(Seq1):
                print ('skipped',Seq0,Seq1,len(Seq0),len(Seq1))
                return(9999999999999999999999999999999999999)				
            else:				
             while c<Len:
                if Seq0[c]=='A' and Seq1[c]=='T': Dif['FP']+=1
                elif Seq0[c]=='T' and Seq1[c]=='A': Dif['FN']+=1 
                elif Seq0[c]=='T' and Seq1[c]=='T': Dif['TP']+=1  
                elif Seq0[c]=='A' and Seq1[c]=='A': Dif['TN']+=1                
                if Seq0[c]!=Seq1[c] and Seq0[c]!='?' and Seq1[c]!='?': Dif['Tot']+=1
                c+=1
             return Dif 
def CountDifNum_RmMiss(Seq0,Seq1):
            Len=len(Seq0)
         		
            Dif=0
            c=0
            if len(Seq0)!=len(Seq1):
                print ('skipped',Seq0,Seq1,len(Seq0),len(Seq1))
                return(9999999999999999999999999999999999999)				
            else:				
             while c<Len:
                if Seq0[c]!=Seq1[c] and Seq0[c]!='?' and Seq1[c]!='?': Dif+=1
                c+=1
             return Dif
def GetIdenSeq(TarSeq,SeqDic):
    Iden=[]
    for i in SeqDic:
        Seq=SeqDic[i]
        Dif=CountDifNum_RmMiss(TarSeq,Seq)	
      	
        if Dif==0: Iden.append(i)
    return Iden	 
def GetBestMatSeq(THVarSeqDic,CellSeq):  

        Dif2THls={}
        for TH in THVarSeqDic:
            C=CountDifNum_RmMiss(CellSeq,THVarSeqDic[TH])
            Dif2THls[C]=Dif2THls.get(C,[])+[TH]
        DifLs=list(Dif2THls.keys())
        DifLs.sort()
     
        BestLs=Dif2THls[DifLs[0]]
        Dic={}
        for B in BestLs:
            Dic[B]=THVarSeqDic[B]
        return BestLs,Dic      
def GetBest(THVarSeqDic,CellSeq,Var):
        print (THVarSeqDic)
        TP2THls={}        
        for TH in THVarSeqDic:
            Dif=CountDifNum_RmMiss_AT(THVarSeqDic[TH],CellSeq)
            TP2THls[Dif[Var]]=TP2THls.get(Dif[Var],[])+[TH]            
        TPLs=list(TP2THls.keys())
        TPLs.sort()
   
        if Var=='TP' or Var=='TN': BestPos=TPLs[-1]  
        else: BestPos=TPLs[0]  
        BestLs=TP2THls[BestPos]
        return TP2THls,list(set(BestLs)) ,BestPos   
      
def SeqAnnoST(THfas,THclo,STcellFas,MutPFile):
    HitCloLs=[]
    Col2Val=ListCol(THclo)
    CloLs=list(set(Col2Val['Node']+['Outgroup']))
    THls,TH2Seq=ReadFasSeq(THfas)
    out=[]
    for Clo in CloLs:
        out.append('>'+Clo+'\n'+TH2Seq['>'+Clo]+'\n')
    THfas1=THfas[:-6]+'_prune.fasta'
    GetOut(THfas1,''.join(out))
    Col2Val=ListCol(MutPFile)
    MutP=list(map(int,Col2Val['label']))
    print (len(MutP)) 
    SeqNum,LabFile,MOAin1=MakeMOAin2(MutP,'',STcellFas)
    SeqAnno(THfas1,STcellFas[:-6]+'MOAIn.fasta','n')
   
def SeqAnno(THFas,OriFas,PosID):
  
    THCellLs,THVarSeqDic=ReadFasSeq(THFas)    
    OriCellLs,OriVarSeqDic=ReadFasSeq(OriFas)
    
    out=['Cell\tNode\tTP\tFP\tFN\tTN\tMissCount']
    for Cell in OriVarSeqDic:
        CellSeq=OriVarSeqDic[Cell]
        BestTHls,BestTHdic=GetBestMatSeq(THVarSeqDic,CellSeq)
        TP2THls,TPBestLs,TPC=GetBest(BestTHdic,CellSeq,'TP')
        FP2THls,FPBestLs,FPC=GetBest(BestTHdic,CellSeq,'FP')
        FN2THls,FNBestLs,FNC=GetBest(BestTHdic,CellSeq,'FN')
        TN2THls,TNBestLs,TNC=GetBest(BestTHdic,CellSeq,'TN')
        BestT=[]
        TH2FP= InvertDic3(FP2THls) 
        for TB in TPBestLs:
            if TB in FPBestLs: BestT.append(TB)
        if BestT==[]: 
            BestT=TPBestLs            
        TH2FN= InvertDic3(FN2THls)     
        FN2bestTH={}
        for BT in BestT:
           FN=TH2FN[BT][0]
           FN2bestTH[FN]=FN2bestTH.get(FN,[])+[BT]
        FNls=list(FN2bestTH.keys())
        FNls.sort()
        BestT1= FN2bestTH[FNls[0]]
        TH2TN= InvertDic3(TN2THls) 
        TNC=TH2TN[BestT1[0]][0]
        MissC=CellSeq.count('?')    
        if len(BestT1) ==1:     
            if PosID!='y':  out.append('\t'.join([Cell,';'.join(BestT1).replace('\"','').replace('>',''),str(TPC),str(TH2FP[BestT1[0]][0]),str(FNls[0]),str(TNC),str(MissC)]))  
            else: 
                 obsFP=TH2FP[BestT1[0]][0]
                 obsFN=FNls[0]
                 ErrorM=1.0*obsFP/(TPC)
                 ErrorW=1.0*obsFN/(TPC)
             
                 out.append('\t'.join([Cell,';'.join(BestT1).replace('\"','').replace('>',''),str(TPC),str(TH2FP[BestT1[0]][0]),str(FNls[0]),str(TNC),str(MissC)]))             
        else:
            if PosID!='y':  out.append('\t'.join([Cell,';'.join(BestT1).replace('\"','').replace('>',''),str(TPC),str(TH2FP[BestT1[0]][0]),str(FNls[0]),str(TNC),str(MissC)])) 
            else: 
                 obsFP=TH2FP[BestT1[0]][0]
                 obsFN=FNls[0]
                 ErrorM=1.0*obsFP/(TPC)
                 ErrorW=1.0*obsFN/(TPC)
  
                 out.append('\t'.join([Cell,';'.join(BestT1).replace('\"','').replace('>',''),str(TPC),str(TH2FP[BestT1[0]][0]),str(FNls[0]),str(TNC),str(MissC)]))      
  
    GetOut(THFas[:-6]+'_CellAnnoAll.txt','\n'.join(out)) 
def SeqAnnoBestMat(THFas,OriFas):  
    THCellLs,THVarSeqDic=ReadFasSeq(THFas)    
    OriCellLs,OriVarSeqDic=ReadFasSeq(OriFas)
    out=['Cell\tNode\tTP\tFP\tFN\tTN\tMissCount']
    for Cell in OriVarSeqDic:
        Dif2THls={}
        for TH in THVarSeqDic:
            C=CountDifNum_RmMiss(OriVarSeqDic[Cell],THVarSeqDic[TH])
            Dif2THls[C]=Dif2THls.get(C,[])+[TH]
        DifLs=list(Dif2THls.keys())
        DifLs.sort()
  
        BestLs=Dif2THls[DifLs[0]]
     
        if len(BestLs)==1: 
            
            if Cell.replace('>','').replace('#','')!='Normal': 
                  out.append(Cell.replace('>','').replace('#','')+'\t'+BestLs[0].replace('>','').replace('#',''))   
    GetOut(THFas[:-6]+'_CellAnnoAll.txt','\n'.join(out))               
def CellAnnotate(GV,THAnno,MutIDFile,OriFas): 
    if OriFas[-3:]=='meg': OriCellLs,OriCell2Seq=ReadMegSeq(OriFas)
    else: OriCellLs,OriCell2Seq=ReadFasSeq(OriFas)
    
    TH2CellLs=ReadTHAnno(THAnno,1)    
    Cell=TH2CellLs[list(TH2CellLs.keys())[0]][0]
    if OriCellLs.count('>'+Cell)==0 and OriCellLs.count('#'+Cell)==0:
        TH2CellLs=ReadTHAnno(THAnno,0) 
    MutPosLs,THVarSeqDic=GV2Fas1(GV,MutIDFile)

    OriVarSeqDic=MakeVarPosSeq(MutPosLs,OriCell2Seq)
    out=['Seq\tTopHapID']
    outG=['#MEGA\n!Title SNVs;\n!Format datatype=dna;\n']
    for Cell in OriVarSeqDic:
        Dif2THls={}
        for TH in THVarSeqDic:
            C=CountDifNum_RmMiss(OriVarSeqDic[Cell],THVarSeqDic[TH])
            Dif2THls[C]=Dif2THls.get(C,[])+[TH]
        DifLs=list(Dif2THls.keys())
        DifLs.sort()
    
        BestLs=Dif2THls[DifLs[0]]
    
        if len(BestLs)==1: 
            
            if Cell.replace('>','').replace('#','')!='Normal': 
                  outG.append('#'+Cell.replace('>','').replace('#','')+'_{'+BestLs[0].replace('>','').replace('#','')+'}'+'\n'+  OriCell2Seq[Cell])
                  out.append(Cell.replace('>','').replace('#','')+'\t'+BestLs[0].replace('>','').replace('#',''))  
      
    Len=len(OriCell2Seq[Cell])       
    GetOut(GV[:-3]+'_CellAnnoAll.txt','\n'.join(out)) 
    outG.append('#Normal_{Normal}\n'+('A'*Len))    
    GetOut(GV[:-3]+'_GroupMeg.meg','\n'.join(outG))

    os.system('megacc -a distance_estimation_between_grp_avg_nucleotide.mao -d '+GV[:-3]+'_GroupMeg.meg -o '+GV[:-3]+'_GroupMegDist.meg')    
    os.system('megacc -a infer_NJ_nucleotide.mao -d '+GV[:-3]+'_GroupMegDist.meg -o '+GV[:-3]+'_GroupMegDist.nwk')

def AdjustNoude2Cout(Dic,PriLs):
    NewDic={}
    for N in Dic:
       
       if str(N).find(';')!=-1:
          Nls=N.split(';')
          Hit=[]
          for N0 in Nls:
              if '\"'+N0+'\"' in PriLs or N0 in PriLs:
                  Hit.append(N0)
          if len(Hit)==1:                  
         
                  NewDic[int(Hit[0])]=NewDic.get(int(Hit[0]),0)+Dic[N]
               
          elif len(Hit)>1:
                  print ('multiple clones were selected',Hit,Nls,PriLs)          
                  open('a','r').readlines()
       else: NewDic[int(N)]=NewDic.get(int(N),0)+Dic[N]           
  
    return NewDic    
def AddBackTopHapMutTree(MOACellAnnoFile,TopHapMutTreeFile,MOAfasta,MOAGV,MinCellC,MinCellC1,FP,FN,SNVc,OriFas,MutP,Last): #MOAupGV_TopHap.gv fasta
    Node2THbackrec={}
    MT=open(TopHapMutTreeFile,'r').readlines()[1:]
    for i in MT:
        i=i.split('\t')
        TH=i[0].strip()
        if TH[:3]=='All':
            N=i[2].split(';')[-1].strip()
            Back=i[3].split(';')
            Rec=i[4].split(';')
            Node2THbackrec[N]={'TH':TH,'R':Rec,'B':Back}
    print (Node2THbackrec)        
    Dec2Anc,NodeMut,Node2In,Edge2In=ReadGVmoa(MOAGV)
    TipLs=GetTipLs(Dec2Anc)
    print ('tipls',TipLs)
    
    df = pd.read_csv(MOACellAnnoFile,sep='\t')
    N2coutDF=df['Node'].value_counts()
    N2cout=N2coutDF.to_dict()

    print (N2coutDF)
    KeepN2NewTHid={}
    KeepN2CellLs={}
    AllN2GoodCellLs={}
    for N in N2cout:
            subdf=df[df['Node']==N]
            subdf.set_index('Cell')
            KeepCellLs=[]
            AllGoodCellLs=[]
            TotFP=0
            TotFN=0
            TotTN=0
            TotTP=0
            for Ind in subdf.index:
               if subdf['TP'][Ind]!=0:
                 obsFP=subdf['FP'][Ind]
                 obsFN=subdf['FN'][Ind]
                 TotError=obsFP+obsFN
                 obsFPR=1.0*obsFP/(obsFP+subdf['TN'][Ind])
                
                 obsFNR=1.0*obsFN/(obsFN+subdf['TP'][Ind])
                 ErrorR=1.0*(obsFN+obsFP)/(obsFN+subdf['TP'][Ind])
                 ErrorM=1.0*obsFP/(subdf['TP'][Ind])
                 ErrorW=1.0*obsFN/(subdf['TP'][Ind])
                 if obsFPR<=FP and obsFNR<=FN: AllGoodCellLs.append(subdf['Cell'][Ind])
                 else:
                     if obsFPR<=FP:
                    
                         p=DoConcTest0(subdf['TP'][Ind],obsFN,FN)#(observed_mut_mut,observed_mut_wild,beta)
                    
                         if p<0.05: print ('FN W count is significantly, large, remove',p,subdf['Cell'][Ind],N)
                         else: AllGoodCellLs.append(subdf['Cell'][Ind])
                     elif obsFNR<=FN:
                       
                         p=DoConcTest0(subdf['TN'][Ind],obsFP,FP)#(observed_mut_mut,observed_mut_wild,beta)
                         if p<0.05: print ('FP M count is significantly, large, remove',p,subdf['Cell'][Ind],N)
                         else: AllGoodCellLs.append(subdf['Cell'][Ind])                         
                     else:
                      
                         p=DoConcTest0(subdf['TP'][Ind],obsFN,FN)#(observed_mut_mut,observed_mut_wild,beta)
                         if p<0.05: print ('FN W count is significantly, large, remove',p,subdf['Cell'][Ind],N)
                         else: #AllGoodCellLs.append(subdf['Cell'][Ind])                         
                            
                             p=DoConcTest0(subdf['TN'][Ind],obsFP,FP)#(observed_mut_mut,observed_mut_wild,beta)
                             if p<0.05: print ('FP M count is significantly, large, remove',p,subdf['Cell'][Ind],N)
                             else: AllGoodCellLs.append(subdf['Cell'][Ind])                           
            KeepCellLs=AllGoodCellLs
            if (len(KeepCellLs)>= MinCellC) :           
             
              NewID='Alladd'
              if N in Node2THbackrec: 
                TH=Node2THbackrec[N]['TH']
                if TH[:3]=='All':
                    NewID=TH
              if str(N).find(';')==-1:      
                  KeepN2NewTHid[str(N)]=NewID+'Node'+str(N) 
                  KeepN2CellLs[str(N)]=KeepCellLs
            if AllGoodCellLs!=[]      :AllN2GoodCellLs[str(N)]=AllGoodCellLs
    for Tip in TipLs:
        Tip=Tip.replace('\"','')
        KeepN2NewTHid[Tip]=KeepN2NewTHid.get(Tip,'AlladdNode'+Tip)

    print ('clone count',len(KeepN2NewTHid),KeepN2NewTHid.keys())
    MissP=GetMissPos(SNVc,MutP)
    print (MissP,len(MissP))
    
    if OriFas[-3:]=='meg': OriCellLs1,OriCell2Seq=ReadMegSeq(OriFas)
    else: OriCellLs1,OriCell2Seq=ReadFasSeq(OriFas)
  
    out=['digraph G {\n']
    for N in Node2In:
        N=N.replace('\"','')        
        out.append(N+' [label=\"'+N+'\\n\"]\n')   
    for E in Edge2In:
        E=E.replace('\"','')
        out.append(E+';\n')   
    out.append('}\n')
    GetOut(MOAGV[:-3]+'_TopHap0.gv',''.join(out))      
    if Last=='n': return []
    
    print (len(KeepN2NewTHid))  
    print ('update fasta')
    NodeLs,Node2Seq=ReadFasSeq(MOAfasta)
    out=[]
    for KeepN in KeepN2NewTHid:
        
        out.append('>'+KeepN2NewTHid[KeepN]+'\n'+Node2Seq['>\"'+KeepN+'\"']+'\n')
    out.append('>Outgroup\n'+('A'*len(Node2Seq['>\"'+KeepN+'\"']))+'\n')
    GetOut(MOAGV[:-3]+'_TopHap.fasta',''.join(out))     #tips not added

    print ('update GV')
    out=['digraph G {\n']
    for N in Node2In:
        N=N.replace('\"','')
        R=[]
        B=[]
        if N in Node2THbackrec:
             R=Node2THbackrec[N]['R']
             B=Node2THbackrec[N]['B']
        Lab=KeepN2NewTHid.get(N,N)
        AllMut=[N]        
        if B!=[] and B!=['']:
            Bin='Back '+'\\nBack '.join(B)
            AllMut.append(Bin)
        if R!=[] and R!=['']:
            Rin='Rec '+'\\nRec '.join(R)
            AllMut.append(Rin)        
        out.append(Lab+' [label=\"'+Lab+'\\n'+'\\n'.join(AllMut)+'\"]\n')
   
    for E in Edge2In:
        Nls=E.replace('\"','').split('->')
        NlsNew=[]
        for N in Nls:
            NlsNew.append(KeepN2NewTHid.get(N,N))
        out.append('->'.join(NlsNew)+';\n')
    out.append('}\n')
    GetOut(MOAGV[:-3]+'_TopHap.gv',''.join(out)) 
    out=['Group\tSeq\tTopHapID\n']
    CellDone=[]
    print ('clone count',len(KeepN2NewTHid))
    for N in KeepN2NewTHid:
        CellLs=KeepN2CellLs.get(N,[])
    
        if CellLs!=[] and CellLs!=['']:
          for Cell in CellLs:
            out.append('NA\t'+Cell.replace('>','')+'\t'+KeepN2NewTHid[N]+'\n')
            CellDone.append(Cell)
        else:
           print ('missing',N,CellLs) 
       
    print ('examine hybrid')
    for N in AllN2GoodCellLs:
        CellLs=AllN2GoodCellLs[N]
        if N in KeepN2NewTHid:

           for Cell in CellLs:
             if Cell not in CellDone:
                out.append('NA\t'+Cell.replace('>','')+'\t'+KeepN2NewTHid[N]+'\n')
                
                CellDone.append(Cell)
        elif N.find(';')!=-1:
           Nls=N.split(';')
      
           TipNs=[]
           for TN in TipLs:
              
               TN=TN.replace('\"','')
               if TN in Nls: TipNs.append(TN)
           
           if len(TipNs)==1:
           
               for Cell in CellLs:
                  if Cell not in OriCell2Seq: Nuc=OriCell2Seq[Cell.replace('>','#')][int(TipNs[0])]
                  else: Nuc=OriCell2Seq[Cell][int(TipNs[0])]
                  if Nuc=='T':
                      if Cell not in CellDone:
                         out.append('NA\t'+Cell.replace('>','')+'\t'+KeepN2NewTHid[TipNs[0]]+'\n')
                         CellDone.append(Cell)
   
           HitNs=[]
           for TN in KeepN2NewTHid:
            
                   TN=TN.replace('\"','')
                   if TN in Nls: HitNs.append(TN)
           if len(HitNs)==1:
                   for Cell in CellLs:
                      if Cell not in CellDone:
                         out.append('NA\t'+Cell.replace('>','')+'\t'+KeepN2NewTHid[HitNs[0]]+'\n')
                        
    GetOut(MOAGV[:-3]+'_TopHap_CellAnno.txt',''.join(out))                 

def SeqrchMissingMut(MissPLs,KeepN2CellLs,OriCell2Seq,FP,FN,Min,TipLs) :
    MissP2AddNode={}
    for MissP in MissPLs:
        
        for N in KeepN2CellLs:
            NucLs={}
            CellLs=KeepN2CellLs[N]
            for Cell in CellLs:
                 Cell=Cell.replace('>','')
                 if '>'+Cell in OriCell2Seq:
                     Seq=OriCell2Seq['>'+Cell][MissP]
                 else: Seq=OriCell2Seq['#'+Cell][MissP]
                 NucLs[Seq]=NucLs.get(Seq,0)+1
            MutC=NucLs.get('T',0)
            WildC=NucLs.get('A',0)
            Tot=MutC+WildC 
            if Tot!=0:
               Mpro=1.0*MutC/Tot
               Wpro=1.0*WildC/Tot 
     
            if MutC>=Min :  MissP2AddNode[MissP]=MissP2AddNode.get(MissP,[])+[N]
  
    MissP2AddNode1={}
    for M in MissP2AddNode:
        if len(MissP2AddNode[M])==1: MissP2AddNode1[M]=MissP2AddNode[M]
    print ('miss add including intermediate',MissP2AddNode1)    
    Node2MissPLs=InvertDic3(MissP2AddNode1) 
    print (Node2MissPLs)
    NewEdge=[]
    NewNode=[]    
    for AncN in Node2MissPLs:
         MissPLs=Node2MissPLs[AncN]  
         if '\"'+AncN+'\"' in TipLs:
            print (AncN,MissPLs)         
            if len(MissPLs)==1 :         
               NewEdge.append(AncN+'->'+str(MissPLs[0]))
               NewNode.append(str(MissPLs[0]))
            else:
                CellLs=KeepN2CellLs[AncN] 
                Pos2CellLs={}
                for Miss1 in MissPLs:
                    for Cell in CellLs:
                        if Cell not in OriCell2Seq:
                             Nuc=OriCell2Seq[Cell.replace('>','#')][Miss1]
                        else:  Nuc=OriCell2Seq[Cell][Miss1]    
                        if Nuc=='T': Pos2CellLs[Miss1]=Pos2CellLs.get(Miss1,[])+[Cell]
                print (Pos2CellLs)
                Cell2Pos=InvertDic3(Pos2CellLs)
                print (Cell2Pos)
                GoodCellC=0
                for Cell in Cell2Pos:
                     if len(Cell2Pos[Cell])==len(MissPLs) :   GoodCellC+=1
                Pos2CellC={}
                for Pos in  Pos2CellLs:
                     Pos2CellC[Pos]=len(Pos2CellLs[Pos]) 
                CellC2Pos=InvertDic1(Pos2CellC)                     
                CellCOr= list(CellC2Pos.keys())
                CellCOr.sort(reverse=True)
                PosOr=[]
                for CellC in CellCOr:
                    PosOr+=CellC2Pos[CellC] 
                print (PosOr)                       
                if GoodCellC>=1:   
                    MissPosC=len(PosOr)
                    MissP1=0
                    while MissP1<MissPosC:
                       if MissP1==0:                    
                          NewEdge.append(AncN+'->'+str(PosOr[MissP1]))
                          NewNode.append(str(PosOr[MissP1]))  
                       else:
                          NewEdge.append(str(PosOr[MissP1-1])+'->'+str(PosOr[MissP1]))
                          NewNode.append(str(PosOr[MissP1]))  
                       MissP1+=1                          
    print (NewEdge,NewNode)
     
    return NewEdge,NewNode
def GetMissPos(SNVc,MutP):
    MissPosLs=[]
    Pos=0
    while Pos<SNVc:
        if Pos not in MutP: MissPosLs.append(Pos)
        Pos+=1 
    return MissPosLs         
def MakeVarPosSeq(PosLs,SeqDic):
    NewDic={}
    for i in SeqDic:
       Seq=SeqDic[i]
       New=''
       for P in PosLs:
          New+=Seq[P]
       NewDic[i]=New
    return (NewDic)
def TestIntermed(ExpFas,OriFas,FP,FN): #ExpFas[:-6] '_FilInter.fasta'    
    BadInterNLs=[]
    CellAnno=ExpFas[:-6]+'_CellAnnoAll.txt'  
    GV=ExpFas[:-6]+'.gv'
    OutFas=ExpFas[:-6]+'_FilInter.fasta'
    Dec2Anc,NodeMut,Node2In,Edge2In	=ReadGV1(GV)
    if OriFas[-3:]=='meg': OcellLs,Ocell2Seq=ReadMegSeq(OriFas)
    else: OcellLs,Ocell2Seq=ReadFasSeq(OriFas)
    EcellLs,Ecell2Seq=ReadFasSeq(ExpFas)
 
    for D in Dec2Anc:
        if D[0]=='A':
            A=Dec2Anc[D]
            if A[0]=='A':
                Pval=DoConcTest(int(A.split('Node')[1]),int(D.split('Node')[1]),Ocell2Seq,FN)
                if Pval<0.01: print ('intermediate detected',A)
                else:
                               
                   BadInterNLs.append(A)
    print ('remove unsupported intermendiate',BadInterNLs)
    out=[]
    for i in EcellLs:
        if i.replace('>','') not in BadInterNLs: out.append(i+'\n'+Ecell2Seq[i]+'\n')
      
    GetOut(OutFas,''.join(out))    
    return BadInterNLs    
  
def DoConcTest(target_variant1,target_variant2,Original_allSequenceList,beta):   #anc dec pos
    result = get_observedFrequency(Original_allSequenceList, target_variant1, target_variant2)  
    observed_mut_mut = len(result[0])
    observed_mut_wild = len(result[1])
    MWpro=1.0*observed_mut_wild/(observed_mut_wild+observed_mut_mut)
    if MWpro<=beta or observed_mut_wild<2: return 1.0
    else:
        result2 = get_expectedFrequency(observed_mut_mut,observed_mut_wild, beta )

        expected_mut_mut = result2[0]
        expected_mut_wild = result2[1]    
        ct = get_contingencyTable(observed_mut_mut, observed_mut_wild, expected_mut_mut, expected_mut_wild)

        oddsr, p = fisher_exact(ct)

        return p  
def DoConcTest0(observed_mut_mut,observed_mut_wild,beta):   #beta=FN

        result2 = get_expectedFrequency(observed_mut_mut,observed_mut_wild, beta )

        expected_mut_mut = result2[0]
        expected_mut_wild = result2[1]    
        ct = get_contingencyTable(observed_mut_mut, observed_mut_wild, expected_mut_mut, expected_mut_wild)

        oddsr, p = fisher_exact(ct)

        return p         
def get_observedFrequency(seq_list, target_v1, target_v2):
    observedMM = []
    observedMW = []	
    for Cell in seq_list:
        sequence=seq_list[Cell]

        if sequence[target_v1] =='T' and sequence[target_v2] == 'A':
            observedMW.append(sequence)

        if sequence[target_v1] =='T' and sequence[target_v2] == 'T':
            observedMM.append(sequence)
    return  observedMM, observedMW

def get_expectedFrequency(observedMM, observedMW, beta):
    total = observedMM + observedMW
    expectedMW = total * beta
    expectedMM = total - expectedMW
    return expectedMM, expectedMW

def get_contingencyTable(oMM, oMW, exMM, exMW):
    table = np.array([[oMM, oMW],[exMM, exMW]])
    return table    
def TestDoublet(CellAnno,ExpFas,OriFas,FPR,FNR): #both Fas is hap
 
    Fol=ExpFas.replace(ExpFas.split(os.sep)[-1],'')
    Cells,Cell2Seq=ReadFasSeq(OriFas)
    ExpCells,ExpCell2Seq=ReadFasSeq(ExpFas)
    UpNode2CellLs={}
    Node2Cell=ReadTHAnno1(CellAnno,1,0)
    for N in Node2Cell:
      if N.find(';')==-1:    
        CellLs=Node2Cell[N]
        ExpSeq=ExpCell2Seq['>'+N]
        GoodCellLs=[]
        DoubPair=[]
        for Cell in CellLs:
               if '>'+Cell not in Cell2Seq: ObsCellSeq=Cell2Seq['#'+Cell]
               else: ObsCellSeq=Cell2Seq['>'+Cell]
               if os.path.exists(Fol+Cell+'.fasta')==True: os.remove(Fol+Cell+'.fasta')
           
               GetFasSeq({Cell:ObsCellSeq}, Fol+Cell+'.fasta')
               IncPos=GetFP(ExpSeq,ObsCellSeq)
               Doub='n'  
                           
               for Pos in IncPos:
                  TarExpSeqDic=GetSeqWithMut(ExpCell2Seq,[Pos])
                  if os.path.exists( Fol+'TestExpSeq.fasta')==True: os.remove( Fol+'TestExpSeq.fasta')
                  GetFasSeq(TarExpSeqDic, Fol+'TestExpSeq.fasta')
              
                  if TarExpSeqDic!={}:  
                    SeqAnno( Fol+'TestExpSeq.fasta',Fol+Cell+'.fasta','n')
                    HitAltSeqLs=open( Fol+'TestExpSeq_CellAnnoAll.txt','r').readlines()[1].split('\t')[1].split(';')
               
                    for HitAlt in HitAltSeqLs:
                       if HitAlt!=N: 
                           HitAltSeq=ExpCell2Seq['>'+HitAlt]
                           print ('\ntest doublet',N,HitAlt,Cell,Pos)
                  
                           LRTpVal=DoLRT(ExpSeq,HitAltSeq,ObsCellSeq,FPR,FNR)
                           print ('p-val',LRTpVal)
                           if float(LRTpVal)<0.01: 
                               Doub='y'
                               DoubPair.append(HitAlt)
                            
               if Doub=='n':GoodCellLs.append(Cell)             
         
        UpNode2CellLs[N]= {'CellLs':GoodCellLs,'pair':list(set(DoubPair))}  
        print (UpNode2CellLs)

    out=['Cell\tNode\n']
    outFas=[]
    for N in UpNode2CellLs:
        CellLs=UpNode2CellLs[N]['CellLs']
        outFas.append('>'+N+'\n'+ExpCell2Seq['>'+N]+'\n')
        if CellLs!=[]:
          for Cell in CellLs:
            out.append('\t'.join([Cell,N])+'\n')
    outFas.append('>Outgroup\n'+('A'*len(ExpCell2Seq['>'+N]))+'\n')
    outFas=list(set(outFas))    
    GetOut(CellAnno[:-4]+'_final.txt',''.join(out))   
    GetOut(CellAnno[:-4]+'_final.fasta',''.join(outFas))  
    
    return (UpNode2CellLs)        
def DoLRT(expected_seq1,expected_seq2,Obs,alpha, beta):
        hybrid = make_hybrid(expected_seq1, expected_seq2)

        l1 = get_ML(Obs, expected_seq1,alpha, beta)
        l2 = get_ML(Obs, hybrid,alpha, beta)
        p = LR_pval(l1,l2)  
        return p
def make_hybrid(exp_1, exp_2):
    hybrid =[]
    for i in range(len(exp_1)):
        if exp_1[i] == 'T' and exp_2[i] != 'T':
            hybrid.append(exp_1[i])
        else:
            hybrid.append(exp_2[i])
    return ''.join(hybrid)        
def get_ML(seq1, seq2,alpha, beta):

    total_prob = []
    c1=0
    c2=0
    for i in range(len(seq1)):
        if seq1[i] == 'A' and seq2[i] == 'A':
            prob = 1 - alpha
        if seq1[i] == 'T' and seq2[i] == 'A':
            prob = alpha
            c1+=1
        if seq1[i] == 'T' and seq2[i] == 'T':
            prob = 1 - beta
        if seq1[i] == 'A' and seq2[i] == 'T':
            prob = beta
            c2+=1
        if seq1[i] == '?' and seq2[i] == 'T':
           prob = beta#1 #- beta
        if seq1[i] == '?' and seq2[i] == 'A':
            prob = alpha#1 #- alpha

        total_prob.append(prob)
    print (c1,c2,alpha,beta)    
    result = np.prod(total_prob)
    return result

def LR_pval(L1,L2):
  print (L1,L2)
  if L1>L2: return 1.0
  else: 
    LR = 2 * abs((np.log(L1) - np.log(L2)))
    p = chi2.sf(LR, 1)
    return ' %.30f' % p                       
def GetSeqWithMut(C2S,PosIncLs):
    NewDic={}
    for C in C2S:
        S=C2S[C]
        Inc='y'
        for P in PosIncLs:
            if S[P]!='T': Inc='n'
        if Inc=='y': NewDic[C]=S 
    return NewDic        
                   
def GetFP(ESeq,OSeq):
    FPs=[]
    Len=len(ESeq)
    c=0
    if len(ESeq)!=len(OSeq): 
        print (len(ESeq),len(OSeq))
        open('a','r').readlines()
    while c<Len:
       if ESeq[c]=='A' and OSeq[c]=='T': FPs.append(c)
       c+=1
    return FPs       
def TopHapAnnotation(TopHap,AllHap):

    Out=TopHap[:-6]+AllHap.split(os.sep)[-1][:-6]+'_anno.txt'
    
    out='Group\tSeq\tTopHapID\n'
    TopLs,TopSeq=ReadFasSeq(TopHap)
    HapLs,HapSeq=ReadFasSeq(AllHap)
    c=0
    for Hap in HapLs:
        HapS=HapSeq[Hap]
        IdenLs=GetIdenSeq(HapS,TopSeq)
  	
        c+=1	
        if len(IdenLs)==1:
     	
            out+='NA\t'+Hap.replace('>','')+'\t'+IdenLs[0].replace('>','')+'\n'
        elif len(IdenLs)>1: pass
	
    GetOut(Out,out)  
def BootstrapSeq_afterBoo(BooC,FillPer,OriTreeFile,OriFas,Rpath,Fol,FasTaLs,HFabsMinMax):
    OutFolMain='Bootstrap' + os.sep#'C:\\Users\\kumarlab\\Desktop\\TopHap\\68k_SNV72\\Boo\\'
    OutGFas=Fol+'OutG.fasta'#'OutGHap.fasta'
    PosDic={} #pos from 1
    Boo2FasLs={}
    go='n'
    BooFolLs=[]
    print (FasTaLs)
    for FasLs in FasTaLs:
        if FasLs==Fol+'Note-Continent.txt'	: Name='MinorCountry_'
        else: Name=''
        FasLs=open(FasLs,'r').readlines()[1:]
        for i in FasLs:
                print ('making bootstrap TopHap seqs',i)
                Fas0=i.strip().split('\t')[0]
                OutFol=OutFolMain+Fas0.replace(' ','_')+os.sep
                if os.path.exists(OutFol)!=True: os.mkdir(OutFol)
                Fas=Fol+Name+Fas0.replace(' ','_')#+'_Hap.fasta'
                StLs,St2Seq=ReadFasSeq(Fas)
                Len=len(StLs)
                Fill=FillPer*Len
                Boo=1

                while Boo<=BooC:
                    Samp=list(np.random.choice(StLs,replace=True,size=Len))
                    Seq2Cou={}
                    for i in Samp:
                       Seq2Cou[St2Seq[i]]=Seq2Cou.get(St2Seq[i],0)+1
                    out=''
                    IDTop=1
                    TopHap2Data=[]					
                    for Seq in Seq2Cou:
                         if Seq2Cou[Seq]>=Fill:
                              out+='>'+str(IDTop)+'\n'+Seq+'\n'
                              TopHap2Data.append(Seq)							  
                              IDTop+=1					
                    CouLs=InvertDic2Ls(Seq2Cou)
          
                    CouLs.sort(reverse=True)
                    print (CouLs)
                    THc=len(TopHap2Data)
                    c=0			
                    while THc<HFabsMinMax[0]:
                         Cou=CouLs[c]
                         for Seq in Seq2Cou:
                          if Seq2Cou[Seq]==Cou:				 
     
                            if Seq not in TopHap2Data:				 
                              out+='>'+str(IDTop)+'\n'+Seq+'\n'
                              TopHap2Data.append(Seq)							  
                              IDTop+=1	
                         THc=len(TopHap2Data)
                         print (Cou,THc)				 
                         c+=1						
                    GetOut(OutFol+'Rep'+str(Boo)+'TopHap.fasta',out)
                    Boo+=1
    print ('make fasta for each boo')
    OriFas_pru=OriFas[:-6]+'_prune.fasta'
    AllRepTrees=''
    Boo=1
    RepFasLs=[]

    while Boo<=BooC:
      print (Boo,'make fasta for each boo')
      TSeq2C={}
      for FasLs in FasTaLs:
        FasLs=open(FasLs,'r').readlines()[1:]
        for i in FasLs:
              
                Fas0=i.strip().split('\t')[0]
                TopHapSeq=OutFolMain+Fas0.replace(' ','_')+os.sep	+'Rep'+str(Boo)+'TopHap.fasta'
                Tls,T2Seq=ReadFasSeq(TopHapSeq)
                for T in T2Seq:
                     Seq=T2Seq[T]
                     TSeq2C[Seq]=TSeq2C.get(Seq,0)+1

      ID=1
      out=''
      for Hap in TSeq2C:
          out+='>A'+str(ID)+'_'+str(TSeq2C[Hap])+'\n'+Hap+'\n'
          ID+=1
      if os.path.exists(OutFolMain+'Res'+str(Boo)+'.fasta') !=True: GetOut(OutFolMain+'Res'+str(Boo)+'.fasta',out)
      RepFasLs.append(OutFolMain+'Res'+str(Boo)+'.fasta')
      Boo+=1
    print ('boo C',len(RepFasLs))
    StLs,St2Seq=ReadFasSeq(OriFas)
    Hap2ID=InvertDic(St2Seq)
    Hap2Count={}
    for St in StLs:

        Hap2Count[St]=0
    OutLs,Out2Seq=ReadFasSeq(OutGFas)
    OutG=''
    for i in OutLs:
        OutG+=i+'\n'+Out2Seq[i]+'\n'

    for RepFas in RepFasLs:
        RStLs,RSt2Seq=ReadFasSeq(RepFas)
        Rhap2ID=InvertDic(RSt2Seq)

        for Rhap in Rhap2ID:
            if Rhap in Hap2ID: Hap2Count[Hap2ID[Rhap]]+=1
    ExtraTip=[]
    out=OutG
    for Hap in Hap2Count:
        if Hap2Count[Hap]!=len(RepFasLs) and OutLs.count(Hap)==0: ExtraTip.append(Hap.replace('>','').replace('__','_'))
        elif OutLs.count(Hap)==0: out+=Hap+'\n'+St2Seq[Hap]+'\n'
    GetOut(OriFas_pru,out)
    print ('extra tip',ExtraTip,len(ExtraTip))
    prune_tree(OriTreeFile,ExtraTip,OutLs[0].replace('>',''))
    for RepFas in RepFasLs:
        RStLs,RSt2Seq=ReadFasSeq(RepFas)
        Rhap2ID=InvertDic(RSt2Seq)
        out=OutG
        ID=RepFas.split(os.sep)[-1].split('_')[0]
        IDC=1
        ExtraLs=[]
        for Rhap in Rhap2ID:
            if Rhap in Hap2ID:
               NewID=Hap2ID[Rhap].replace('>','').replace('__','_')
               if ExtraTip.count(NewID)!=0: ExtraLs.append(NewID)
            else:
                NewID=ID+'_'+str(IDC)
                ExtraLs.append(NewID)
                IDC+=1

            out+='>'+NewID+'\n'+Rhap+'\n'
        if os.path.exists(RepFas[:-6]+'_Rename.fasta')!=True:
              GetOut(RepFas[:-6]+'_Rename.fasta',out)

        if os.path.exists(RepFas[:-6]+'_Rename.nwk')!=True:
          os.system('megacc -a infer_MP_nucleotide.mao -d ' + RepFas[:-6]+'_Rename.fasta' +' -o ' + RepFas[:-6]+'_Rename.nwk')
          AncFileLs=glob.glob(RepFas[:-6]+'_Rename_ancestral_states_*.txt')
          for AncF in AncFileLs:
               os.remove(AncF)

        print (ID,len(ExtraLs),len(RStLs),len(RStLs)-len(ExtraLs))
        if os.path.exists(RepFas[:-6]+'_Rename.nwk')==True: prune_tree(RepFas[:-6]+'_Rename.nwk',ExtraLs,OutLs[0].replace('>',''))
        PruTreeLs=glob.glob(RepFas[:-6]+'_Rename_*_prune.nwk')
        for PruTre in PruTreeLs:
            TrLines=open(PruTre,'r').readlines()
            for Tr in TrLines:
                AllRepTrees+=Tr

    GetOut('AllResampTrees.nwk',AllRepTrees)
    OriPruneLs=glob.glob(OriTreeFile[:-4]+'_*_prune.nwk')
    print ('target tres',OriPruneLs)
    trees = list(Phylo.parse('AllResampTrees.nwk', "newick"))
    c=1
    for OriTreFile in OriPruneLs:
        print ('score tree',OriTreFile,ID)
        BooSum(OriTreFile,str(c),Rpath)

        c+=1
        os.remove(OriTreFile)


def FillMissBase_consen(Fas,GoodFilLs,SupCut):

    Out=Fas[:-6]+'_FillAllCell.fasta'
    CellLs,Cell2Seq=ReadFasSeq(Fas)
    print ('tot cell c',len(CellLs))
    CellWithMiss=[]
    for Cell in CellLs:
        Seq=Cell2Seq[Cell]
        if Seq.find('?')!=-1 or Seq.find('-')!=-1:
            CellWithMiss.append(Cell)
    print ('cell with missing base',len(CellWithMiss))
    out=''
    Rec=0
    CellC=0
    TotCell=len(CellLs)
    Unique=[]
    print ('filling missing data...')
    for Cell in CellLs:

        CellC+=1	
        Seq=Cell2Seq[Cell]
        if CellWithMiss.count(Cell)==0:
             out+=Cell+'\n'+Cell2Seq[Cell]+'\n'
        else:

           MatchHap2C={}
           for Cell1 in CellLs:
               if Cell!=Cell1:
                    Seq1=Cell2Seq[Cell1]		   
                    DiffNum=CountDifNum_excMiss(Seq,Seq1)
                    if DiffNum==0:
                        MatchHap2C[Seq1]=Cell2Seq.get(Seq1,0)+1

           if MatchHap2C=={}:
               Unique.append(Cell+'\n'+Seq)	   
           Pos2Fill={}
           c=0
           Len=len(Seq)
           Filled=''	   
           UnFillC=0
           FillC=0	   
           while c<Len:
               if Seq[c]=='?' or Seq[c]=='-':
                   Fill={}
                   TotMatch=0			   
                   for Mhap in MatchHap2C:
                       Nuc=Mhap[c]
                       if Nuc!='?' and Nuc!='-':
                             Fill[Nuc]=Fill.get(Nuc,0)+MatchHap2C[Mhap]
                             TotMatch+=MatchHap2C[Mhap]						 

                   Cou2NucLs=InvertDic1(Fill)
                   Cou=list(Cou2NucLs.keys())			   
                   Cou.sort()
                   if len(Cou)==0: 
                        Filled+='?'			   
                        UnFillC+=1  
			
                   else:					
                    BestLs=Cou2NucLs[Cou[-1]]
                    Sup=1.0*Cou[-1]/TotMatch				
			
                    if len(BestLs)==1 and GoodFilLs.count(list(Fill.keys())[0])!=0 and Sup>=SupCut:
		   
                        Filled+=BestLs[0]
                        FillC+=1					
                    else: 
                        Filled+='?'			   
                        UnFillC+=1 	
					
               else: Filled+=Seq[c]						
               c+=1
           if UnFillC==0: 
               out+=Cell+'\n'+Filled+'\n'
               Rec+=1
           elif FillC!=0:
               out+=Cell+'\n'+Filled+'\n'	   
               Rec+=1	
           else: 
               out+=Cell+'\n'+Filled+'\n'	   
               Rec+=1		   
		   
    GetOut(Out,out)   
    GetOut('UniqueSeq.fasta','\n'.join(Unique))    	   
def AnnotateOriID(THanno,cellIDanno):  

    Out=THanno[:-4]+'1.txt'
    cellIDanno=open(cellIDanno,'r').readlines()[1:]
    THID2OriID={}
    for i in cellIDanno:
        i=i.split('\t')
        THID2OriID[i[0]]=i[3]

    THanno=open(THanno,'r').readlines()
    out=THanno[0]
    THanno=THanno[1:]
    for i in THanno:
        i=i.split('\t')
        TH=i[1]
        if TH.split('_')[0] in THID2OriID:	
           Ori=THID2OriID[TH.split('_')[0]]	
           out+=TH+'\t'+Ori+'\t'+i[2].replace('__','_')	

    GetOut(Out,out) 
def TopHap_DrawMutTree2(Fas,PosID,Tree,OutGroup):

    PosIDorder=open(PosID,'r').readlines()
    GetOut('outgroup.txt',OutGroup+'=outgroup')
    
    OutF=Fas[:-6]+'_MutTree.txt'
    
    root_tree(Tree,OutGroup)	
    
    if os.path.exists(Tree[:-4]+'_Anc.csv')!=True:
        os.system('megacc -a ancestral_seqs_ML_nucleotideNucJK.mao -d '+Fas+' -t '+Tree[:-4]+'_rooted.nwk -o '+Tree[:-4]+'_Anc.csv -g outgroup.txt')

    MakeAncSeq(Tree[:-4]+'_Anc.csv',OutGroup)
    AncLs,Anc2Seq=ReadMegSeq(Tree[:-4]+'_AncWithAnc.meg')
    D2A,A2D,N2C,C2N=ReadNodeMap(Tree[:-4]+'_Anc_nodeMap.txt')
    
    Root=''
    for A in A2D:
        if A not in D2A: Root=A
	
    OutC=N2C[OutGroup]
    OutGNodeLs=[OutC]

    Bra2Info={}
    Old2Young={} #zero branch
    for N in N2C:
        if N[:5]!='Node_':
           Anc=N2C[N]
           if Anc==Root: pass
           else:           	   
            Bot=N	   
            BotSeq=Anc2Seq['#Node_'+N2C[N]+'CellID'+N]	   
            while Anc in D2A:
               Old=Anc	   
               Anc=D2A[Anc]
               if Anc==Root:
                 Seq=Anc2Seq['#Node_'+OutC+'CellID'+C2N[OutC]]		   
               else:		   
                 Seq=Anc2Seq['#Node_'+Anc+'CellID'+C2N[Anc]]
               DifC=CountDifNum(BotSeq,Seq)
               if DifC>0:	
                   if Bot not in Bra2Info:	
                   
                      ForLs,BackLs=ClassifyDif1(Seq,BotSeq,PosIDorder)				  
                      Bra2Info[Bot]=[C2N[Anc],ForLs,BackLs]  			   
                   Bot=C2N[Anc]
                   BotSeq=Seq
               else: Old2Young[Anc]=Old 			   
 
    Mut2C={}
    
    Bra2Info1={}
    for Bra in Bra2Info:
        InfoLs=Bra2Info[Bra][1]
        for i in InfoLs:
            Mut2C[i]=Mut2C.get(i,0)+1
    	
        Bra1=N2C[Bra]		
        if Bra1 in Old2Young:
            Young=Old2Young[Bra1]
            while Young in Old2Young:
                Young=Old2Young[Young]
            Bra1=Young
 
      
        Anc1=N2C[Bra2Info[Bra][0]]
        if Anc1 in Old2Young:
           Young=Old2Young[Anc1]
           while Young in Old2Young:
                Young=Old2Young[Young]
           Anc1=Young
    
        Bra2Info1[C2N[Bra1]]=[C2N[Anc1]]+Bra2Info[Bra][1:]	
    
    In2DecLs={}
    Mut2C={}
    for Bra in Bra2Info1:
        InfoLs=Bra2Info1[Bra][1]
        for i in InfoLs:
            Mut2C[i]=Mut2C.get(i,0)+1
        Info=Bra2Info1[Bra]
        In='\t'.join(map(str,Info))
        In2DecLs[In]=In2DecLs.get(In,[])+[Bra]
    AncID2IDup={}	
    RmLs=[]
    for In in In2DecLs:
         DecLs=In2DecLs[In]
         if len(DecLs)>1:
         
            IDup=DecLs[0]
            for i in DecLs:
               if i[:5]!='Node_': IDup=i	
            for i in DecLs:
               if i!=IDup: 
                 AncID2IDup[i]=IDup
                 RmLs.append(i)			 
   
    out='Branch\tAncestor\tMuttation\tLoss of mutation\tParallel mutations\tCell count\n'
    for Bra in Bra2Info1:
        Info=Bra2Info1[Bra]
        For=Info[1]
        if Info[2]!=[]: 
         
           open('a','r').readlines()	   

        if RmLs.count(Bra)==0:	
           	
            out+=Bra+'\t'+AncID2IDup.get(Info[0],Info[0])+'\t'+';'.join(map(str,For))+'\t'+';'.join(map(str,Info[2]))+'\t\n'

    GetOut(OutF,out)

    MutTree2dot1(OutF)
    Dec2Anc,NodeMut,Node2In,Edge2In=ReadGV(OutF[:-4]+'.gv')
    RootLs=[]
    ChangeN2Root=[]
    Anc2Dec=GetAnc2Dec(Dec2Anc)
    for A in Anc2Dec:
        if A not in Dec2Anc:
         if A not in NodeMut: ChangeN2Root.append(A)
         else:	 
            MutLs=NodeMut[A]
            Mc=0
            for M in MutLs: 
                if M.strip()!='': Mc+=1
            if Mc==0: ChangeN2Root.append(A)			
            else: RootLs.append(A)
    print ('root',RootLs,ChangeN2Root)
    out='digraph D {\n'
    for N in Node2In:
        if N in ChangeN2Root: pass
     
        else: out+=Node2In[N]
    out+='root [label=\"root\"]\n'
 
    for E in Edge2In:
        Eori=E.split('->')[0].strip()
        if Eori in 	ChangeN2Root: out+='root'+Edge2In[E][len(Eori):]+'\n'
     
        else: out+=Edge2In[E]+'\n'	
     
    for R in RootLs:
          out+='root->'+R+';\n'	
    out+='}'
    GetOut(OutF[:-4]+'.gv',out) 
def GetBestNuc(Dic):
    TMP={}
    for Node in Dic:
        Nuc2Prob=Dic[Node]
        BestP=0
        BestNuc='?'
        for Nuc in Nuc2Prob:
            P=Nuc2Prob[Nuc]
            if P>BestP:
                BestNuc=Nuc
                BestP=P
        TMP[Node]=BestNuc
    return TMP  
def GetCOIID(Mut):	
        Pos=int(Mut[1:][:-1])#+1
        M=Mut[-1]
        MutID=str(Pos)+M
        return MutID
def GetAncLs(MutLs):
    UpID=[]
    for i in MutLs:
      if i!='NA' and i.strip()!='':	
       UpID.append(GetCOIID(i))
    return UpID	
def GetAncMut(Anc,NodeMut,Dec2Anc):
    Ls=[]
    while Anc in Dec2Anc:
       Anc=Dec2Anc[Anc]
       if Anc in NodeMut: Ls+=NodeMut[Anc]
    return Ls	
def GetDecMut(Dec0,NodeMut,Dec2Anc):
   Anc2Dec={}
   for Dec in Dec2Anc:
      Anc=Dec2Anc[Dec]
      Anc2Dec[Anc]=Anc2Dec.get(Anc,[])+[Dec]	  
   Decs=Anc2Dec.get(Dec0,[])
   Add='y'
   while Add=='y':
      AddLs=[]
      Len1=len(Decs)	  
      for i in Decs:
         AddLs+=Anc2Dec.get(i,[])
      AddLs=list(set(AddLs))		 
      Decs+=AddLs	
      Decs=list(set(Decs))
      Len2=len(Decs)	  
      if Len1==Len2: Add='n'
   Ls=[]
   for i in Decs:
       if i in NodeMut: Ls+=NodeMut[i]   
	
   return (Ls)
   
def GetRec(MutDic):
    Mut2C={}
    for N in MutDic:
        MutLs=MutDic[N]
        for M in MutLs:
            Mut2C[M]=Mut2C.get(M,0)+1
    C2MutLs=InvertDic1(Mut2C)
    RecLs=[]
    for C in C2MutLs:
       if C!=1 :
        Cand=C2MutLs[C]
        for M in Cand:
          if M.strip()!='':		
             RecLs.append(M)	
    return RecLs,Mut2C

	
def SumCOI(COI,GV,COIcut):

    Dec2Anc2COI=ReadCOI(COI)
    Pat2COILs=ListCol(COI) #child/parent 0T

    
    Dec2Anc,NodeMut,Node2In,Edge2In=ReadGV(GV)
    Anc2Dec=GetAnc2Dec(Dec2Anc)


    InLs=['Type\tTH1\tTH2\tRelationship\tMut1\tMut2\tCOI\tSD\tPairCount\tPvalTtes']
    Node2RecBack={}
    NodeLs1=[]
    NodeLs2=[]    
    for Node in NodeMut:
             MutLs=NodeMut[Node]

             for M in MutLs:
               if M.strip()!='' and M.find('?')==-1:
                if M[0]!='A': Node2RecBack[Node]=Node2RecBack.get(Node,[])+[M]
                else:
                 
                   for M1 in MutLs:
                     if M1.strip()!='' and M1.find('?')==-1:
                       if M!=M1 and M1[0]=='A':
                            COI=Dec2Anc2COI[GetCOIID(M)][GetCOIID(M1)]
                            InLs.append('\t'.join(['For',Node,Node,'DecAnc(Conc)',M,M1,str(COI),'NA','NA','NA']))
                            NodeLs1.append(Node+'\t'+Node)
                            NodeLs2.append(Node)                            
                     
                   Anode=Node 
                                      
                   while Anode in Dec2Anc:
                        Anode=Dec2Anc[Anode]
                        AMutLs=NodeMut[Anode] 
                        for AM in AMutLs:
                        
                            if AM.strip()!='' and AM.find('?')==-1:
                            
                              if AM[0]=='A':
                                COI=Dec2Anc2COI[GetCOIID(M)][GetCOIID(AM)]
                                InLs.append('\t'.join(['For',Node,Anode,'DecAnc',M,AM,str(COI),'NA','NA','NA']))
                                NodeLs1.append(Node+'\t'+Anode)
                                NodeLs2.append(Anode)                                 
 
    for Node in Node2RecBack:
        RecBackLs=Node2RecBack[Node]
        DecNodeLs=GetDecCla(Node,Dec2Anc)
        for RB in RecBackLs:
            Type=''
            if RB[0]=='R': Type='Rec'
            elif  RB[0]=='T': Type='Back'
            for Dnode in DecNodeLs:
                        DMutLs=NodeMut[Dnode] 
                        for DM in DMutLs:
                          
                            if DM.strip()!='' and RB.find('?')==-1:
                            
                              if DM[0]=='A':
                                COI=Dec2Anc2COI[GetCOIID(DM)][GetCOIID(RB.replace('Rec',''))]
                                InLs.append('\t'.join([Type,Dnode,Node,'DecAnc',DM,RB,str(COI),'NA','NA','NA'])) 
                                NodeLs1.append(Dnode+'\t'+Node)
                                NodeLs2.append(Node)                                 
    GetOut(GV.replace(GV.split(os.sep)[-1],'')+'COIsummary.txt','\n'.join(InLs)) 
    import pandas
    import scipy
    data=pandas.read_csv(GV.replace(GV.split(os.sep)[-1],'')+'COIsummary.txt', sep='\t')
    AvePat=data.groupby(['Type','TH2','Mut2']).agg({'COI':['mean','std','count']})
    CombPat,OrderPat=SumPanda(AvePat)
    Mut2Pval={}
       
    for TypeTHmu in CombPat:
            sub3=data[(data['Mut2']==TypeTHmu[2]) & (data['TH2']==TypeTHmu[1])]
  
            if len(sub3['COI'])<5:
              
                    if  sub3['COI'].mean()> COIcut: pvalue=0
                    else: pvalue=1 
            else:  
                Ttesres=scipy.stats.ttest_1samp(sub3['COI'],COIcut, alternative='greater')  
                pvalue=Ttesres.pvalue                
            Item=[TypeTHmu[0],'-',TypeTHmu[1],'-','-',TypeTHmu[2],CombPat[TypeTHmu][0],CombPat[TypeTHmu][1],CombPat[TypeTHmu][2],pvalue]
            InLs.append('\t'.join(map(str,Item)))
            if TypeTHmu[2] in Mut2Pval: Mut2Pval[TypeTHmu[2]][TypeTHmu[1]]={'DecSupp':pvalue}
            
            else: Mut2Pval[TypeTHmu[2]]={TypeTHmu[1]:{'DecSupp':pvalue}}
    for Node in NodeMut: 
      if Node in Dec2Anc:
        DirectAnc=Dec2Anc[Node]
        if NodeLs1.count(Node+'\t'+Node)!=0 :
         
           sub=pandas.read_csv(GV.replace(GV.split(os.sep)[-1],'')+'COIsummary.txt', sep='\t')
        
           sub1=sub[sub['TH1']==Node]
            
           if NodeLs1.count(Node+'\t'+DirectAnc)!=0 :  
                        
              sub2=sub1[(sub1['TH2']==Node) | (sub1['TH2']==DirectAnc)]#sub1.loc[[DirectAnc,Node]]
          
           else:  sub2=sub1[sub1['TH2']==Node]#sub1.loc[Node]  
           AveChi=sub2.groupby(['Type','Mut1']).agg({'COI':['mean','std','count']})
      
           CombChi,OrderChi=SumPanda(AveChi)
                 
           for TypeTHmu in CombChi:  
             if  TypeTHmu[0]=='For': 
                sub3=sub2[(sub2['Mut1']==TypeTHmu[1]) & (sub2['Type']=='For')]
                
                if len(sub3['COI'])<5:
                    if sub3['COI'].mean()> COIcut: pvalue=0
                    else: pvalue=1   
                else:  
                    Ttesres=scipy.stats.ttest_1samp(sub3['COI'],COIcut, alternative='greater')
                    pvalue=Ttesres.pvalue                    
                Item=[TypeTHmu[0],Node,'-','-',TypeTHmu[1],'-',CombChi[TypeTHmu][0],CombChi[TypeTHmu][1],CombChi[TypeTHmu][2],pvalue]
                InLs.append('\t'.join(map(str,(Item))))
                if TypeTHmu[1] not in Mut2Pval: Mut2Pval[TypeTHmu[1]]={}
                if Node not in Mut2Pval[TypeTHmu[1]]:  Mut2Pval[TypeTHmu[1]][Node]={}
                Mut2Pval[TypeTHmu[1]][Node]['AncSupp']=pvalue
            
    GetOut(GV.replace(GV.split(os.sep)[-1],'')+'COIsummary1.txt','\n'.join(InLs))                      
    FilterMutByCOI(Mut2Pval,GV)   
def SumPanda(Ave):
    Ave1=Ave.to_dict()
    Comb={}
    for MeanStd in Ave1:
      
        In=Ave1[MeanStd]
     
        for i in In:
            Comb[i]=Comb.get(i,[])+[In[i]] 
    return (Comb,list(Ave1))
def FilterMutByCOI(Mut2Pval,GV):
    Out=GV[:-3]+'_ave3.gv'
    Dec2Anc,NodeMut,Node2In,Edge2In=ReadGV(GV)
    Anc2Dec=GetAnc2Dec(Dec2Anc)
    out='digraph D {\n'
    for Node in Node2In:
        MutLs=NodeMut[Node]
        MutIn=[]
        for M in MutLs:
            if M not in Mut2Pval:
                 print (M,'no COI support, so removed')
            elif Node not in Mut2Pval[M]:  print (M,Node,'no COI support, so removed')    
            else:      
                COIpDic=Mut2Pval[M][Node]
                Good='n'
                Bad='n'
                for COIp in COIpDic:
                    Val=COIpDic[COIp]
                    if Val!='nan':
                        Val=float(Val)
                        if Val<=0.01: Good='y'
                        else: Bad='y'
                if Good=='y' and Bad=='n': MutIn.append(M)
                else: print (M,'low COI support, so removed') 
        Good='n'
        for M in MutIn:
           if M[0]=='A': Good='y'
        if Good=='n' and Node not in Anc2Dec:
           print (Node,MutIn,'all mutations are recurrent/back and it is Tip, so removed')
           MutIn=[]           
        out+= Node+' [label=\"'+Node+'\\nMut:'+';'.join(MutIn)+'\\n\"]\n'
    for E in Edge2In:        
        out+=E+'\n'
    out+='}\n'
    GetOut(Out,out)    
    
def MapCOI3(COI,GV):
 
    Dec2Anc2COI=ReadCOI(COI)
    Pat2COILs=ListCol(COI) #child/parent 0T
    TrunkLs=[]
    TrunkLsA=[]
    for Pat in Pat2COILs:
        if Pat!='child\parent':
           COILs=Pat2COILs[Pat]
           Trunk='y'
           for i in COILs:
              if float(i)<0.6: Trunk='n'
           if Trunk=='y': 
              TrunkLs.append(Pat)
              TrunkLsA.append('A'+Pat)		  
    print ('trunk mut',TrunkLs)	   
    
    
    Dec2Anc,NodeMut,Node2In,Edge2In=ReadGV(GV)
    Anc2Dec=GetAnc2Dec(Dec2Anc)
    RootDecLs=Anc2Dec['root']
    if len(RootDecLs)==1: Root=RootDecLs[0]
    else: Root='root'
    print ('root node',Root)
    NodeMut[Root]=list(set(NodeMut[Root]+TrunkLsA))
    Node2In[Root]=Root+' [label=\"'+Root+'\\nMut:'+';'.join(NodeMut[Root])+'\\n\"]\n'
    for Node in NodeMut:
        if Node!=Root:#'root':
             MutLs=NodeMut[Node]
             New=[]
             for M in MutLs:
                if 	TrunkLsA.count(M)==0: New.append(M)
                else: print ('prune from NodeMut',Node,M)
             NodeMut[Node]=New			
    
       
    RecLs,Mut2C=GetRec(NodeMut)
    print ('recurrent mut',RecLs)

       
    out='digraph D {\n'
    for Node in Node2In:
        In=Node2In[Node]
        for Rec in RecLs:
           In=In.replace(Rec,'Rec'+Rec)	
        if Node!=Root:#'root':	   
         for Tru in TrunkLs:
          if  In.find('A'+Tru)!=-1:
              print ('pruned trunk mut',Node,Tru)
              In=In.replace('A'+Tru,'')			  
        out+=In#Node2In[Node]
    AllCOI=[]
	
    for Edge0 in Edge2In:
        Edge=Edge0.split('->')
        Anc=Edge[0].strip()
        Dec=Edge[1].strip()

        AncMutLs=NodeMut.get(Anc,['NA'])
        DecMutLs=NodeMut[Dec]
        AllAncMutLs=GetAncMut(Anc,NodeMut,Dec2Anc)	
        AllDecMutLs=GetDecMut(Dec,NodeMut,Dec2Anc)	

        AncLs=GetAncLs(AncMutLs)+GetAncLs(DecMutLs)+GetAncLs(AllAncMutLs)+GetAncLs(AllDecMutLs)#[AMid,DMid]

        In=[]
	
        for DM in DecMutLs:
          if DM.strip()!='':# and DM[-1]=='T':	
           DMid=GetCOIID(DM)
           if TrunkLs.count(DMid)==0:
	   
            if DMid not in Dec2Anc2COI	or RecLs.count('A'+DMid)!=0:
                 In.append(DM+':NA')
                 Anc2Coi={}			 
		 
            else: 		
               Anc2Coi=Dec2Anc2COI[DMid]	
    	   
               for AM in AncMutLs:
                if AM=='NA': In.append(DM+':NA')
                else:
                 if AM.strip()!='':			
                  AMid=GetCOIID(AM)
		  
                  COI=Anc2Coi.get(AMid,'NA')
                  if COI=='NA':
                       In.append(DM+':'+AM+'=NA')
                  else:
                       AllCOI.append(COI)			  

                       In.append(DM+':'+AM+'='+str(COI))
				   
            if RecLs.count('A'+DMid)==0:				   
             Best=[]
             BestCOI=0
	 
             for A in Anc2Coi:
                 B='A'+A		 
                 if AncLs.count(A)==0:# and RecLs.count(B)==0:
                      if B in Mut2C:				   
                           C=Anc2Coi[A]
                           if C>BestCOI:
                               BestCOI=C	
                               Best=[A]
                           elif C==BestCOI: Best.append(A)
             if BestCOI==0: In.append('/'+DM+':NA')	
             else: In.append('/'+DM+':'+','.join(Best)+'='+str(BestCOI))
		
        out+=Edge0+' [label=\"  '+'\\n'.join(In)+'\"];\n'
    out+='}\n'

    GetOut(GV[:-3]+'_COI3.gv',out)    
def AveCOI3(GV):
    AltFil='n'
    AltCut=0.7
    COICut=0.4
  
    Out=GV[:-3]+'_ave3.gv'
    Dec2Anc,NodeMut,Node2In,Edge2In	=ReadGV(GV)
    out='digraph D {'
    for Node in Node2In:
        out+=Node2In[Node]+'\n'
    AllAve=[]
    AllAlt=[]
    PruneLs=[]
    AltH2C={}	
    for Edge in Edge2In:
       In=Edge2In[Edge]
       PreChi=''  
       if In.find('\"')==-1: out+=In+'\n'
       else:
        InLs=In.split('\"')[1].split('\\n')
        Chi2Pat2COI={}
        Chi2Pat2Alt={}	
        for i in InLs:
          if i.strip()!='':	
            i=i.strip().split(':')
         	
            Chi=i[0]
          		
            PatCOI=i[1].split('=')
            Pat=PatCOI[0]
        	
            if Pat!='NA':
             COI=PatCOI[1]#.split('/')[0]
             if COI!='NA': 
               COI=float(COI)		
               if Chi[0]!='/':		   
    			
                Chi2Pat2COI[Chi]=Chi2Pat2COI.get(Chi,[])+[COI]
                PreChi=Chi			
               else: 
                 if COI>AltCut: 
                     AltH2C[Pat]=AltH2C.get(Pat,0)+1
                  				 
                 else: Chi2Pat2Alt[Chi[1:]]=Chi2Pat2Alt.get(Chi[1:],[])+[Pat+':'+str(COI)]#.append(COI)#[Chi]=Chi2Pat2COI.get(Chi,[])+[COI]
    		 
    		   
        ChiCOI=[]
      
        if Chi2Pat2COI=={}:
           out+=Edge+' [label=\"NA/'#
           if Chi2Pat2Alt=={}: out+='NA\"];\n'
           else: 
              AltLs=[]	   
              SupportLs=[]	
              PatAltIn=[]			  
              for Chi in Chi2Pat2Alt:	
                 if len(Chi2Pat2Alt[Chi])!=1: open('a','r').readlines()
                 else:
    		 
                    for PatCOI in Chi2Pat2Alt[Chi]:
                   				
                      if Chi[-1]=='T':				
                        AltLs.append(float(PatCOI.split(':')[-1]))
                   
                        PatAltIn.append(Chi+':'+PatCOI.split(':')[0]+':'+str(round(float(PatCOI.split(':')[-1]),2))+'\\n')
              out+='NA/'.join(PatAltIn)+'\"];\n'
              				
              if len(AltLs)!=0:                    					
                AltAve=sum(AltLs)/len(AltLs)
            	
                AllAlt.append(AltAve)
              if len(SupportLs)!=0:           			
                SuppAve=sum(SupportLs)/len(SupportLs)
        
                AllAve.append(SuppAve)			  
        else:
         out+=Edge+' [label=\"'	
         In=[]	
         AltLs=[]	
    	 
         for Chi in	Chi2Pat2COI:
    
                COIls=Chi2Pat2COI[Chi]	
                Ave=sum(COIls)/len(COIls)
                if Ave<COICut:
    
                    print ('prune',Edge,Chi,Pat,COI)			
                    PruneLs.append([Chi,Edge.split('->')[1]])
                else:				
                 PatAlt=Chi2Pat2Alt.get(Chi,[])
                 if PatAlt==[]:In.append(Chi+':'+str(round(Ave,2))+'/NA')
                 elif len(PatAlt)!=1: open('a','r').readlines()
                 else:
                   PatAlt=PatAlt[0].split(':')
                   AltVal=float(PatAlt[1])			   
              		   
                   In.append(Chi+':'+str(round(Ave,2))+'/'+PatAlt[0]+':'+str(round(AltVal,2)))	
                   AltLs.append(float(PatAlt[1]))			   
                 ChiCOI.append(Ave)
    	
         if len(ChiCOI)>0:	
          Ave=sum(ChiCOI)/len(ChiCOI)
          AllAve.append(Ave)	 
          if AltLs!=[]: 	 
             AllAlt.append(sum(AltLs)/len(AltLs))	 
      	 
         out+='\\n'.join(In)+'\"];\n'#str(Ave)+'/'+str(AveAlt)+'\"];\n'
    out+='}\n'	
    GetOut(Out,out)
    Ave=sum(AllAve)/len(AllAve)
    if len(AllAlt)==0: AltAve='NA'
    else:
       AltAve=sum(AllAlt)/len(AllAlt)
    print (Ave,AltAve) 
    Med=numpy.median(AllAve)
    MedA=numpy.median(AllAlt)
    print (Med,MedA)
    
    print ('prune dure to COI<',COICut,PruneLs)
 
    out=''
    AltHLs=[]
    if AltFil=='y':
      for i in AltH2C:
        AltHLs+=i.split(',')
        out+='AltCOI>'+str(AltCut)+': '+i+'\t'+str(AltH2C[i])+'\n'
    print (out)
    AltHLs=list(set(AltHLs))
    print (AltHLs)	
    GetOut(GV.replace(GV.split(os.sep)[-1],'prune.txt'),'\n'.join(str(PruneLs))+'\n'+out)
    AllRm=PruneLs#+AltHLs
    print ('remove ls: ',AllRm)

    GVori=open(Out,'r').readlines()
    out=''
    for i in GVori:
        for RNode in AllRm:
         R=RNode[0]
         Node=RNode[1]	
         if i.split(' ')[0]==Node:	 

            if R[0]!='A' and R[0]!='T': 
                if R[-1]=='T': R='A'+R
                elif R[-1]=='A': R='T'+R
                else: open('a','r').readlines()			
            i=i.replace(R+';','').replace(R,'')
            i=i.replace(';;','')
        out+=i		
    GetOut(GV.replace(GV.split('_')[-1],'')+'prune.gv',out)	
    GetOut(Out,out)		

    GetOut(GV[:-3]+'_AveCOI3.txt','COI: '+str(Ave)+'\nAlt: '+str(AltAve)+'\nCOIMed: '+str(Med)+'\nAltMed: '+str(MedA))   
    os.system('dot -Tpng '+Out+' -o '+Out[:-3]+'.png')     			
def GetComAnc(NodeLs,Dec2Anc):
    Node2C={}
    for N in NodeLs:
        Node2C[N]=Node2C.get(N,0)+1
        while N in Dec2Anc:
            N=Dec2Anc[N]
            Node2C[N]=Node2C.get(N,0)+1
    Com=[]
    for N in Node2C:			
        if Node2C[N]>1: Com.append(N)
    return Com
def getCOI(TarM,AncMLs,Dec2Anc2COI):
    AncCOI=[]
    for A in AncMLs:
      if A!='':	
       if TarM.replace('Rec','')[1:]!='' and A.replace('Rec','')[1:]!='':	
        if TarM.replace('Rec','')[1:] in Dec2Anc2COI:
           COI0=Dec2Anc2COI[TarM.replace('Rec','')[1:]]
           if A.replace('Rec','')[1:] in COI0:           
               COI=Dec2Anc2COI[TarM.replace('Rec','')[1:]][A.replace('Rec','')[1:]]
               AncCOI.append(A+'='+str(COI))
    return AncCOI
def ReadTHAnno(File,CellPos):	
    File=open(File,'r').readlines()[1:]
    TH2CellLs={}
    for i in File:
        i=i.split('\t')
        TH=i[2].strip()
        TH2CellLs[TH]=TH2CellLs.get(TH,[])+[i[CellPos]]
    return TH2CellLs 
def ReadTHAnno1(File,THPos,CellPos):	
    File=open(File,'r').readlines()[1:]
    TH2CellLs={}
    for i in File:
        i=i.split('\t')
        TH=i[THPos].strip().replace('>','').replace('#','')
        
        TH2CellLs[TH]=TH2CellLs.get(TH,[])+[i[CellPos].replace('>','').replace('#','')]
    return TH2CellLs       
def CompFPFN(CellFas,THAnno):
    THfas=THAnno.replace(THAnno.split(os.sep)[-1],'')+'TopHap.fasta'
    THls,TH2Seq=ReadFasSeq(THfas)
    if CellFas[-4:]=='.meg':CellLs,Cell2Seq=ReadMegSeq(CellFas)
    else: CellLs,Cell2Seq=ReadFasSeq(CellFas)
    TH2CellLs=ReadTHAnno(THAnno,1)
    Cell=TH2CellLs[list(TH2CellLs.keys())[0]][0]
    if CellLs.count('>'+Cell)==0 and CellLs.count('#'+Cell)==0:
        TH2CellLs=ReadTHAnno(THAnno,0)   
    Cout={}
    TH2count={}
    for TH in TH2CellLs:
       TH2count[TH]={}
       THseq=TH2Seq['>'+TH]
       Len=len(THseq)
       Cells= TH2CellLs[TH]
       for Cell in Cells:
          CellSeq=Cell2Seq['>'+Cell]
          if len(CellSeq)!=Len: open('a','r').readlines()
          c=0
          while c<Len: 
             Pair=THseq[c]+CellSeq[c]
             Cout[Pair]=Cout.get(Pair,0)+1 
             TH2count[TH][Pair]=TH2count[TH].get(Pair,0)+1 
             c+=1 
    FP=1.0*Cout.get('AT',0)/(Cout.get('AT',0)+Cout.get('AA',0)+Cout.get('A?',0))
    FN=1.0*Cout.get('TA',0)/(Cout.get('TA',0)+Cout.get('TT',0)+Cout.get('T?',0))       
    return FP,FN ,TH2count   
def GetInter(NodeLs,Dec2Anc):
    ComAnc=GetComAnc(NodeLs,Dec2Anc)

    Inter=[]
    Node2Inter={}
    for N in NodeLs:
        A=Dec2Anc[N]
        Node2Inter[N]=[]
        while ComAnc.count(A)==0: 
           Inter.append(A)
           Node2Inter[N].append(A)#=Node2Inter.get(N,[])+[A]
           if A not in Dec2Anc: A=ComAnc[0]
           else: A=Dec2Anc[A]
    return (Inter,Node2Inter)
def GetNuc(Rec,CellLs,Cell2SeqOri):
    Pos=int(Rec[1:][:-1])
    print (Pos,Rec,CellLs)
    Nucs=[]
    print (len(Cell2SeqOri[list(Cell2SeqOri.keys())[0]]))
    for Cell in CellLs:
      if Cell.replace('>','').replace('#','')!='Normal' and Cell.replace('>','').replace('#','')!='Outgroup':     
        if '>'+Cell in Cell2SeqOri: Nucs.append(Cell2SeqOri['>'+Cell][Pos])
        elif '#'+Cell in Cell2SeqOri: Nucs.append(Cell2SeqOri['#'+Cell][Pos])
        else:
            print (Cell,list(Cell2SeqOri.keys())[:3],len(list(Cell2SeqOri.keys())))
            open('a','r').readlines()        
    return (Nucs)    
def GetDecTH(B,N,Dec2Anc,NodeMut):
    print (B)
    Pos=int(B.strip()[1:][:-1])
    print (B,Pos)
    Anc2Dec=InvertDic1(Dec2Anc) 
    DecLs=Anc2Dec.get(N,[])
    print (N,DecLs)
    GoodDec=[]
    while DecLs!=[]:
        NextDecLs=[]
        for Dec in DecLs:
            MutLs=NodeMut[Dec]
            Good='y'
            for M in MutLs:
                if M=='A'+str(Pos)+'T': Good='n'
            if Good=='y': 
                 print ('h',GoodDec)
                 GoodDec.append(Dec)
                 if Dec in Anc2Dec: NextDecLs+=Anc2Dec[Dec]
            print (Dec,NextDecLs,MutLs,Good)
        DecLs=NextDecLs 
    print ('Res',B,N,GoodDec)        
    return GoodDec    
def GetDecCla(N,Dec2Anc):

    Anc2Dec=InvertDic1(Dec2Anc) 
    DecLs=Anc2Dec.get(N,[])
   
    AllDec=[N]
    if N in Anc2Dec:
        DecLs=[N]
        while DecLs!=[]:
           NextDecLs=[]
           for Dec in DecLs:

                 if Dec in Anc2Dec: NextDecLs+=Anc2Dec[Dec]
        
           DecLs=NextDecLs
           AllDec+=NextDecLs            
       
    return AllDec    
def GetInferBack(B,N,Dec2Anc,NodeMut):
    Pos=int(B.replace('Rec','')[1:][:-1])
    print (B,Pos)
    Inter=[]
    A=Dec2Anc[N]
    Mls=NodeMut[A]
    if Mls.count('A'+str(Pos)+'T')==0:
        Inter.append(A)
        Go='y'
        while Go=='y':
             if A not in Dec2Anc: Go='n'
             else:
                 A=Dec2Anc[A]
                 if NodeMut[A].count('A'+str(Pos)+'T')!=0: Go='n'
                 else: Inter.append(A)
    ForCla=GetDecCla(A,Dec2Anc)
    BackCla=GetDecCla(N,Dec2Anc)
    print (N,B,ForCla,BackCla)
    Inter1=[]
    for F in ForCla: 
        if BackCla.count(F)==0: Inter1.append(F)
    print (Inter1)    
  
    return Inter1,A  
def GetAncTH(N,Dec2Anc):
    Als=[]
    while N in Dec2Anc:
       N=Dec2Anc[N]
       Als.append(N)
    return Als   

def AnaTipRec(GV,COI,HapFas,FP,FN):
    TopHapAnno=GV.replace(GV.split(os.sep)[-1],'')+'TopHapNA__2022_1_FillAllCell_anno1.txt'
    TopHap2CellLs=ReadTHAnno(TopHapAnno,1)
    print (HapFas)
    if HapFas[-4:]=='.meg': CellLs,Cell2SeqOri=ReadMegSeq(HapFas)
    else: CellLs,Cell2SeqOri=ReadFasSeq(HapFas)
    Cell=TopHap2CellLs[list(TopHap2CellLs.keys())[0]][0]
    if CellLs.count('>'+Cell)==0 and CellLs.count('#'+Cell)==0:
        TopHap2CellLs=ReadTHAnno(TopHapAnno,0)      
    
    Out=GV.replace(GV.split(os.sep)[-1],'')+'MutTree_Final.gv'

    Dec2Anc,NodeMut,Node2In,Edge2In	=ReadGV1(GV)
    Dec2Anc2COI=ReadCOI(COI)
    
    Rec2NodeLs={}
    for N in NodeMut:
        MutLs=NodeMut[N]
     
        for M in MutLs:
            if M.find('Rec')!=-1 : Rec2NodeLs[M.replace('Rec','').strip()]=Rec2NodeLs.get(M.replace('Rec','').strip(),[])+[N]
  
    Node2RecCOIin={}	
  #  print ('Rec',Rec2NodeLs)
   
    Remove='n'
    Rec2KeepN={}
    for Rec in Rec2NodeLs:
           NodeLs=Rec2NodeLs[Rec]
           Rec=Rec
           Rec2KeepN[Rec]=NodeLs      
           AllDecs1=[]
           AllTip='y'
           FirstSplit='n'
           for N in NodeLs:                            
             DecTHLs=GetDecTH(Rec.replace('Rec','').strip(),N,Dec2Anc,NodeMut) #error do not use
             if len(DecTHLs)>0: AllTip='n'
           if AllTip=='y': 
            #  print ('\n\ndec should have mutation',Rec,N,NodeLs) 
            

              AllDecs1=NodeLs           
              for DecTH in AllDecs1:

                 NucLs=GetNuc(Rec,TopHap2CellLs[DecTH],Cell2SeqOri)
              #   print (DecTH,NucLs)
                 if NucLs!=[]:        
                   MutPro=1.0*NucLs.count('A')/(NucLs.count('T')+NucLs.count('A'))
                 #  print (MutPro,NucLs.count('T'),NucLs.count('A'),FN)    
                   if MutPro<FN : pass# Keep.append(DecTH)
                   else:     
                      Rec2KeepN[Rec]=[]
                      Remove='y'  

                 else: 
                     Rec2KeepN[Rec]=[]  
                     Remove='y' 
                    
  #  print (Rec2KeepN,Remove)        
 
         
    RecNodeLs=InvertDic3(Rec2NodeLs)
    RecKeepNodeLs=	InvertDic3(Rec2KeepN)
    RecNodeLs=list(RecNodeLs.keys())
 
    out='digraph D {\n'
    for N in Node2In:

           In=Node2In[N]
           MutLsO=NodeMut[N]
           KeepRec=RecKeepNodeLs.get(N,[])	   
        
           NewIn=[]
           for OM in MutLsO:
              if OM.find('Back')!=-1: NewIn.append(OM.replace('Back','').strip())
              else:     
               if OM.strip()!='':              
                if OM[:3]!='Rec': NewIn.append(OM)

                elif KeepRec.count(OM[3:])!=0 or KeepRec.count(OM[3:].strip())!=0: NewIn.append('Rec'+OM[3:].strip())
                else: 
                        print ('rec mutation removed',N,OM)		
                                             
         				
           out+=N+'  [label=\"'+N+'\\nMut:'+';'.join(NewIn)+'\\n\"]\n'	   
    for E in Edge2In:
        out+=E+'\n'
    out+='}'
    GetOut(Out,out)  
  #  print (Out)
    return Remove

    
def AnaRecBack(GV,COI,HapFas,FP,FN):
    TopHapAnno=GV.replace(GV.split(os.sep)[-1],'')+'TopHapNA__2022_1_FillAllCell_anno1.txt'
    TopHap2CellLs=ReadTHAnno(TopHapAnno,1)
  #  print (HapFas)
    if HapFas[-4:]=='.meg': CellLs,Cell2SeqOri=ReadMegSeq(HapFas)
    else: CellLs,Cell2SeqOri=ReadFasSeq(HapFas)
    Cell=TopHap2CellLs[list(TopHap2CellLs.keys())[0]][0]
    if CellLs.count('>'+Cell)==0 and CellLs.count('#'+Cell)==0:
        TopHap2CellLs=ReadTHAnno(TopHapAnno,0)  

    Out=GV[:-3]+'_rec.gv'
  
    Dec2Anc,NodeMut,Node2In,Edge2In	=ReadGV(GV)
    Dec2Anc2COI=ReadCOI(COI)
    
    Rec2NodeLs={}
    for N in NodeMut:
        MutLs=NodeMut[N]
        for M in MutLs:
            if M.find('Rec')!=-1 and M[-1]=='T': Rec2NodeLs[M.replace('Rec','')]=Rec2NodeLs.get(M.replace('Rec',''),[])+[N]
    	
    Node2RecCOIin={}	
 #   print (Rec2NodeLs)
    
    Rec2KeepN={}
    Rec2RmN={}
    for Rec in Rec2NodeLs:
      
      if Rec!='':
       NodeLs=list(set(Rec2NodeLs[Rec]))
       if len(NodeLs)==1:
         Rec2KeepN[Rec]=Rec2KeepN.get(Rec,[])+[NodeLs[0]]   
       else:
      	
       	InterThLs,N2InterLs=GetInter(NodeLs,Dec2Anc) #Inter should not have mutations
     #   print (Rec,InterThLs,N2InterLs,NodeLs)
        NoInterNode=[]
        Go='y'
        for TarN in N2InterLs:
          #  print ('\n',TarN,Rec)
            InterLs=N2InterLs[TarN]
            THIc=0
            NucLs=[]      
            for InterTH in InterLs:
                   if InterTH in TopHap2CellLs:
                     NucLs+=GetNuc(Rec,TopHap2CellLs[InterTH],Cell2SeqOri)
                    # print ('inter',InterTH,NucLs)
                   
         #   print (TarN,InterLs,NucLs,NucLs.count('T'),NucLs.count('A'))  
            THIc=len(NucLs)

            if NucLs!=[] and ((NucLs.count('T')+NucLs.count('A'))>0):        
                MutPro=1.0*NucLs.count('T')/(NucLs.count('T')+NucLs.count('A'))
              #  print ('Intermediate found',MutPro,NucLs.count('T'),NucLs.count('A'),FP)    
                if MutPro<FP or NucLs.count('T')==1:  Rec2KeepN[Rec]=Rec2KeepN.get(Rec,[])+[TarN]
                else: Go='n'
                
            else: Rec2KeepN[Rec]=Rec2KeepN.get(Rec,[])+[TarN] ##descendant should have mutations if there is no inter  
           

         #   print (Rec2KeepN)         
            Others=[] #others should not have mutations
          #  print ('others should not have mutations',Others)            
            AllDecs=[]                   
            for N in NodeLs:  
              DecTHLs=GetDecTH(Rec.replace('Rec',''),N,Dec2Anc,NodeMut)
              
              AllDecs+=DecTHLs+[N]            
            for D in Dec2Anc:
                 if AllDecs.count(D)==0 and  InterThLs.count(D)==0 : Others.append(D)

            NucLs=[]        
            for OTH in Others:
               if OTH in TopHap2CellLs:
                 NucLs+=GetNuc(Rec,TopHap2CellLs[OTH],Cell2SeqOri)
              #   print (OTH,NucLs)
                 
          #  print (NucLs,NucLs.count('T'),NucLs.count('A'))  
            if (NucLs.count('T')+NucLs.count('A'))!=0:        
                MutPro=1.0*NucLs.count('T')/(NucLs.count('T')+NucLs.count('A'))
              #  print (MutPro,NucLs.count('T'),NucLs.count('A'))    
                if MutPro<FP or NucLs.count('T')==1:   pass
                else: Go='n'#Rec2KeepN[Rec]=[]                
                    
        if Go=='n': Rec2KeepN[Rec]=[] 
        else: #both tips tend to have false recurrent mutations
           AllDecs1=[]
           AllTip='y'
           FirstSplit='n'
           for N in NodeLs:   
             AncTHLs=GetAncTH(N,Dec2Anc)                          
             DecTHLs=GetDecTH(Rec.replace('Rec',''),N,Dec2Anc,NodeMut)
             if len(DecTHLs)>0: AllTip='n'
             
             AllDecs1+=DecTHLs+[N]
             Cont='y' #first split
             if len(AncTHLs)>2 or len(DecTHLs)>0: Cont='n'
             elif len(AncTHLs)==2:
                 
               #  print (NodeMut[AncTHLs[0]],NodeMut[AncTHLs[1]],AncTHLs)
                 if len(NodeMut[AncTHLs[0]])>0 and len(NodeMut[AncTHLs[1]])>0: Cont='n'
             if Cont=='y':               
                 NucLs=GetNuc(Rec,TopHap2CellLs[N],Cell2SeqOri)
                 if NucLs!=[]:        
                   MutPro=1.0*NucLs.count('A')/(NucLs.count('T')+NucLs.count('A'))
                  # print (MutPro,NucLs.count('T'),NucLs.count('A'),FN)    
                   if MutPro<FN :  pass
                   else:     Rec2KeepN[Rec]=[]#Rec2KeepN.get(Rec,[])+[TarN]#ec2KeepN.get(Rec,[])+[NodeLs[0]]
                 else: Rec2KeepN[Rec]=[]               
             
             
             
           if AllTip=='y': 
            #  print ('\n\ndec should have mutation',Rec,N,DecTHLs,NodeLs) 
            
              AllDecs1=list(set(AllDecs1))
             # print (AllDecs1)
            
              AllDecs1=list(set(AllDecs1))           
              for DecTH in AllDecs1:

                 NucLs=GetNuc(Rec,TopHap2CellLs[DecTH],Cell2SeqOri)
 
                 if NucLs!=[]:        
                   MutPro=1.0*NucLs.count('A')/(NucLs.count('T')+NucLs.count('A'))
                   print (MutPro,NucLs.count('T'),NucLs.count('A'),FN)    
                   if MutPro<FN :  pass#Keep.append(DecTH)
                   else:     Rec2KeepN[Rec]=[]#Rec2KeepN.get(Rec,[])+[TarN]#ec2KeepN.get(Rec,[])+[NodeLs[0]]
                 else: Rec2KeepN[Rec]=[]
                     
            
 #   print (Rec2KeepN)        

         
    RecNodeLs=InvertDic3(Rec2NodeLs)
    RecKeepNodeLs=	InvertDic3(Rec2KeepN)
    RecNodeLs=list(RecNodeLs.keys())
 
    out='digraph D {\n'
    for N in Node2In:
        if N not in RecNodeLs: out+=Node2In[N]
        else:
           In=Node2In[N]
           MutLsO=NodeMut[N]
           KeepRec=RecKeepNodeLs.get(N,[])	   
        
           NewIn=[]
           for OM in MutLsO:
                
                if OM[:3]!='Rec': NewIn.append(OM)
                else:
                    if KeepRec.count(OM[3:])!=0: NewIn.append(OM[3:])
                    else: 
                        print ('rec mutation removed',N,OM)		
                                             
         				
           out+=N+'  [label=\"'+N+'\\nMut:'+';'.join(NewIn)+'\\n\"]\n'	   
    for E in Edge2In:
        out+=E+'\n'
    out+='}'
    GetOut(Out,out)  

    GV=Out
    Out=GV[:-3]+'_back.gv'
 
    Dec2Anc,NodeMut,Node2In,Edge2In	=ReadGV(GV)
 
    
    BackM2NodeLs={}
    for N in NodeMut:
        MutLs=NodeMut[N]
        for M in MutLs:
          if M.strip()!='':	
           if M[-1]=='A': BackM2NodeLs[M]=BackM2NodeLs.get(M,[])+[N]
   # print ('\n\n\n\nback list',BackM2NodeLs)
    	 
    	  
    Node2RmBack={}
    for B in BackM2NodeLs:
        NodeLs=BackM2NodeLs[B]
        for N in NodeLs:
            DecTHls=GetDecTH(B.replace('Rec',''),N,Dec2Anc,NodeMut) #descendant should not be mutated
            DecTHls.append(N)
          #  print ('decendant',B,N,DecTHls)
            NucLs=[]
            for TH in DecTHls:
               if TH in TopHap2CellLs:
                 NucLs+=GetNuc(B.replace('Rec',''),TopHap2CellLs[TH],Cell2SeqOri)
              #   print (TH,NucLs)
                 
          #  print (TH,NucLs,NucLs.count('T'),NucLs.count('A'))  
            if NucLs!=[] and (NucLs.count('T')+NucLs.count('A'))>0:        
                MutPro=1.0*NucLs.count('T')/(NucLs.count('T')+NucLs.count('A'))
               # print (MutPro,NucLs.count('T'),NucLs.count('A'))    
                if MutPro>FP :  Node2RmBack[N]=Node2RmBack.get(N,[])+[B]	 ############remove or NucLs.count('T')==1########
            InterTH,ForNode=GetInferBack(B,N,Dec2Anc,NodeMut)     #inter should be mutated. If not remove both back and forward  
           # print ('inter',B,N,InterTH)
            NucLs=[]
            for TH in InterTH:
               if TH in TopHap2CellLs:
                 NucLs+=GetNuc(B.replace('Rec',''),TopHap2CellLs[TH],Cell2SeqOri)
                # print (TH,NucLs)
                 
           # print (TH,NucLs,NucLs.count('T'),NucLs.count('A'))  
            if NucLs!=[] and (NucLs.count('T')+NucLs.count('A'))>0:        
                MutPro=1.0*NucLs.count('T')/(NucLs.count('T')+NucLs.count('A'))
              #  print (MutPro,NucLs.count('T'),NucLs.count('A'))    
                if MutPro<0.5 and MutPro<FP:  
                    Node2RmBack[N]=Node2RmBack.get(N,[])+[B]	###remove forward also  
                    Node2RmBack[ForNode]=Node2RmBack.get(ForNode,[])+['A'+B[1:][:-1]+'T']                    
            
          # print ('remove back mut',Node2RmBack)
           
 
 #   print ('remove bad backward',Node2RmBack)
    
    out='digraph D {\n'
    for N in Node2In:
        if N not in Node2RmBack: out+=Node2In[N]
        else:
           In=Node2In[N]
           MutLsO=NodeMut[N]
           RmLs=Node2RmBack[N]
           NewIn=[]
           for OM in MutLsO:
                if RmLs.count(OM)==0: NewIn.append(OM)
                else:
                    print ('back mutation removed',N,OM)								
           out+=N+'  [label=\"'+N+'\\nMut:'+';'.join(NewIn)+'\\n\"]\n'	   
    for E in Edge2In:
        out+=E+'\n'
    out+='}'
    GetOut(Out,out)	 
   
def AnaRec1(GV,COI):
   
    Out=GV[:-3]+'_rec.gv'
   
    Dec2Anc,NodeMut,Node2In,Edge2In	=ReadGV(GV)
    Dec2Anc2COI=ReadCOI(COI)
    
    Rec2NodeLs={}
    for N in NodeMut:
        MutLs=NodeMut[N]
        for M in MutLs:
            if M.find('Rec')!=-1: Rec2NodeLs[M.replace('Rec','')]=Rec2NodeLs.get(M.replace('Rec',''),[])+[N]
    	
    Node2RecCOIin={}	
  #  print (Rec2NodeLs)
    
    Rec2KeepN={}
    for Rec in Rec2NodeLs:
      if Rec!='':
       NodeLs=list(set(Rec2NodeLs[Rec]))
       if len(NodeLs)==1:
         Rec2KeepN[Rec]=Rec2KeepN.get(Rec,[])+[NodeLs[0]]   
       else:
        ComAnc=GetComAnc(NodeLs,Dec2Anc)	
       	
        
        for N in NodeLs:
           RecCOIin=[]	
           Anc=N#Dec2Anc[N]
           
           while ComAnc.count(Anc)==0:	   
               COILs=getCOI(Rec,NodeMut[Anc],Dec2Anc2COI)		
               for COI in COILs:
                   if COI.strip()!='':
                      RecCOIin.append(float(COI.split('=')[-1]))			   
               Anc=Dec2Anc[Anc]	
                    
           if len(RecCOIin)<3: 
              P=0	
              Ave=999		  
           else:	    
            Ave=sum(RecCOIin)/len(RecCOIin)
            if Ave>0.1:	   
               A=stats.ttest_1samp(RecCOIin,0.1)
               P=A.pvalue
            else: P=999
           if P<0.05:
               Rec2KeepN[Rec]=Rec2KeepN.get(Rec,[])+[N]
         
    RecNodeLs=InvertDic3(Rec2NodeLs)
    RecKeepNodeLs=	InvertDic3(Rec2KeepN)
    RecNodeLs=list(RecNodeLs.keys())
    
    out='digraph D {\n'
    for N in Node2In:
        if N not in RecNodeLs: out+=Node2In[N]
        else:
           In=Node2In[N]
           MutLsO=NodeMut[N]
           KeepRec=RecKeepNodeLs.get(N,[])	   
        
           NewIn=[]
           for OM in MutLsO:
                if OM[:3]!='Rec': NewIn.append(OM)
                else:
                    if KeepRec.count(OM[3:])!=0: NewIn.append(OM[3:])
                    else: print ('rec mutation removed',N,OM)				
         				
           out+=N+'  [label=\"'+N+'\\nMut:'+';'.join(NewIn)+'\\n\"]\n'	   
    for E in Edge2In:
        out+=E+'\n'
    out+='}'
    GetOut(Out,out)
def TestBackRec(BackRecGV,FinalGV,CellAnno,OriFas,FP,FN):
    OutGV=FinalGV[:-3]+'_backrec.gv'
    rbDec2Anc,rbNodeMut,rbNode2In,rbEdge2In	=ReadGV1(BackRecGV)
  
    fiDec2Anc,fiNodeMut0,fiNode2In0,fiEdge2In	=ReadGV1(FinalGV)
    fiNodeMut={}
    for i in fiNodeMut0:
        In=fiNodeMut0[i]
        New=[]
        for ii in In:
            if ii[0]!='R' and ii[0]!='B' and ii[-1]!='A' and ii.find('?')==-1: New.append(ii)
        fiNodeMut[i]=New
        
    fiNode2In={}
    for NodeID in fiNode2In0:
        MutIn=fiNodeMut[NodeID] 
               
        fiNode2In[NodeID]=NodeID+' [label=\"'+NodeID+'\\n'+'\\n'.join(MutIn)+'\"]\n'        
     
    Pos2NodeID={}
    for fD in fiDec2Anc:
        ID=int(fD.split('Node')[-1])
        Pos2NodeID[ID]=fD
    NodeID2Pos=InvertDic(Pos2NodeID)    
    
    TH2CellLs=ReadTHAnno1(CellAnno,1,0)
    
    if OriFas[-3:]=='meg': CellLs,Cell2Seq=ReadMegSeq(OriFas)
    else: CellLs,Cell2Seq=ReadFasSeq(OriFas)
    RecBacSup=['\t'.join(['Type','Mutation postion (from 0)','NodeID','Mutant Cell Count','Wild Cell Count','Mutant cell proportion','Mutant cell proportion'])]
    for N in rbNodeMut:
       MutLs=rbNodeMut[N]
       RecBaDic=GetRekBack(MutLs)
      
       for RB in RecBaDic:
          if RB!='F':
           MposLs=RecBaDic[RB]
           if len(MposLs)>0:
               ConFLs=RecBaDic['F']
               if len(ConFLs)>0:
                   AllDec=[]
                
                   LenLs=[]
                   Len2F={}
                   for F in ConFLs:
                      Fid=Pos2NodeID[F]
                      DecLs=GetDecCla(Fid,fiDec2Anc)
                      LenLs.append(len(DecLs))
                      Len2F[len(DecLs)]=Len2F.get(len(DecLs),[])+[F]
                   
                      AllDec+=DecLs
                   LenLs.sort()
                 
                   AllDec=list(set(AllDec))
                   MostAnc=Len2F[LenLs[-1]][0]
                   MostDec=Len2F[LenLs[0]]
                 
                   if len(AllDec)<LenLs[-1]: pass#print ('mutation order changed, so recurrent/back mutation cannot be added',MposLs)
                   else: 
                       DecCellLs=[]
                       for TH in AllDec:
                           if TH in TH2CellLs: DecCellLs+=TH2CellLs[TH]
                      
                       for M in MposLs:
                           Go='y'
                           if RB=='B': #does not need this. tree was locked
                               
                               if Pos2NodeID[M] in AllDec: 
                                   
                                      Go='n'
                               else:
                                   DecN=AllDec[0]
                                   Go='n'
                                   while DecN in  fiDec2Anc:
                                        if int(DecN.split('Node')[-1])==M: 
                                             Go='y'
                                            
                                             break 
                                        DecN=fiDec2Anc[DecN]  
                                   if Go=='n': pass#print ('forward mutation not found, so remove',M)                                        
                           if Go=='y':
                               CountDic=CountMut(DecCellLs,Cell2Seq,M) 
                               Tot=CountDic.get('T',0)+CountDic.get('A',0)
                         #      print ('\n',N,RB,M,CountDic)                                
                               if Tot>0:
                                  Mpro=1.0*CountDic.get('T',0)/Tot
                                  Wpro=1.0*CountDic.get('A',0)/Tot
                                  RecBacSup.append('\t'.join([RB,str(M),N,str(CountDic.get('T',0)),str(CountDic.get('A',0)),str(Mpro),str(Wpro)]))
                                  if RB=='B':
                                       if Wpro>FN:
                                           NodeID=Pos2NodeID[MostAnc]

                                           MutIn=fiNodeMut[NodeID]
                                         #  print (MutIn)
                                           MutIn.append('Back '+str(M))
                                           fiNodeMut[NodeID]=MutIn
                                         #  print (NodeID)                                           
                                           fiNode2In[NodeID]=NodeID+' [label=\"'+NodeID+'\\n'+'\\n'.join(MutIn)+'\"]\n'
                                         #  print ('back added',M,NodeID,fiNode2In[NodeID])
                                  elif RB=='R':
                                     #  print (Mpro,FP, Wpro,FN)
                                       if Mpro>FP:
                                           NodeID=Pos2NodeID[MostAnc]

                                           MutIn=fiNodeMut[NodeID]
                                        
                                           MutIn.append('Rec '+str(M))
                                           fiNodeMut[NodeID]=MutIn
                                                                                
                                           fiNode2In[NodeID]=NodeID+' [label=\"'+NodeID+'\\n'+'\\n'.join(MutIn)+'\"]\n'
                                        #   print ('rec added',M,NodeID,fiNode2In[NodeID])                                           

    out=['digraph G {\n']
    for N in  fiNode2In:
        out.append(fiNode2In[N])
    for E in fiEdge2In:
        out.append(E+';\n')
    out.append('}')
    GetOut(OutGV,''.join(out))  
    GetOut(OutGV.replace(OutGV.split(os.sep)[-1],'')+'RecBackSupp.txt','\n'.join(RecBacSup))

def TestBackRec1(BackRecGV,FinalGV,CellAnno,OriFas,FP,FN):
    OutGV=FinalGV[:-3]+'_backrec.gv'
    rbDec2Anc,rbNodeMut,rbNode2In,rbEdge2In	=ReadGV1(BackRecGV)
    rbMut2Node=InvertDic3(rbNodeMut)
#    print (rbMut2Node)

    fiDec2Anc,fiNodeMut0,fiNode2In0,fiEdge2In	=ReadGV1(FinalGV)
    fiNodeMut={}
    for i in fiNodeMut0:
        In=fiNodeMut0[i]
        New=[]
        for ii in In:
            if ii[0]!='R' and ii[0]!='B' and ii[-1]!='A' and ii.find('?')==-1: New.append(ii)
        fiNodeMut[i]=New
 
    fiNode2In={}
    for NodeID in fiNode2In0:
        MutIn=fiNodeMut[NodeID] 
              
        fiNode2In[NodeID]=NodeID+' [label=\"'+NodeID+'\\n'+'\\n'.join(MutIn)+'\"]\n'        
   
    Pos2NodeID={}
    for fD in fiDec2Anc:
        ID=int(fD.split('Node')[-1])
        Pos2NodeID[ID]=fD
    NodeID2Pos=InvertDic(Pos2NodeID)    
    
    TH2CellLs=ReadTHAnno1(CellAnno,1,0)
    
    if OriFas[-3:]=='meg': CellLs,Cell2Seq=ReadMegSeq(OriFas)
    else: CellLs,Cell2Seq=ReadFasSeq(OriFas)
    RecBacSup=['\t'.join(['Type','Mutation postion (from 0)','NodeID','Mutant Cell Count','Wild Cell Count','Mutant cell proportion','Mutant cell proportion'])]
    BM2For={}
    for N in rbNodeMut:
       MutLs=rbNodeMut[N]
       RecBaDic=GetRekBack(MutLs)
      
       for RB in RecBaDic:
          if RB!='F':
           MposLs=RecBaDic[RB]
           if len(MposLs)>0:
               ConFLs=RecBaDic['F']
               if len(ConFLs)>0:
                   AllDec=[]
                 
                   LenLs=[]
                   Len2F={}
                   F2DecLs={}
                   for F in ConFLs:
                     if F in Pos2NodeID:
                      Fid=Pos2NodeID[F]
                      DecLs=GetDecCla(Fid,fiDec2Anc)
                      LenLs.append(len(DecLs))
                      Len2F[len(DecLs)]=Len2F.get(len(DecLs),[])+[F]
                      F2DecLs[F]=DecLs
                   
                      AllDec+=DecLs
                   LenLs.sort()
                
                   AllDec=list(set(AllDec))
                   if len(LenLs)>0:
                     MostAnc=Len2F[LenLs[-1]][0]
                     MostDec=Len2F[LenLs[0]][0]
               
                     if len(AllDec)<LenLs[-1]: pass#print ('mutation order changed, so recurrent/back mutation cannot be added',MposLs)
                     else: 
                       DecCellLs=[]
                       for TH in F2DecLs[MostDec]:
                           if TH in TH2CellLs: DecCellLs+=TH2CellLs[TH]
                    
                       for M in MposLs:
                           Go='y'
                           if RB=='B': 
                            
                               Ancb=fiDec2Anc.get(Pos2NodeID[MostAnc],'')
                               AllAnc=[] #from dec
                               if Ancb!='':
                                   if Ancb!='root': AllAnc.append(NodeID2Pos[Ancb])
                                   while Ancb in fiDec2Anc:
                                         Ancb=fiDec2Anc[Ancb]
                                         if Ancb!='root': AllAnc.append(NodeID2Pos[Ancb])
                               
                               FirstConcFor=GetFirstConcMut(rbNodeMut[rbMut2Node['A'+str(M)+'T'][0]],AllAnc)    
                               if FirstConcFor!='': BM2For['A'+str(M)+'T']=FirstConcFor 
                               else: Go='n'
                                                                 
                           if Go=='y':
                               CountDic=CountMut(DecCellLs,Cell2Seq,M) 
                               Tot=CountDic.get('T',0)+CountDic.get('A',0)
                            #   print ('\n',N,RB,M,CountDic)                                
                               if Tot>0:
                                  Mpro=1.0*CountDic.get('T',0)/Tot
                                  Wpro=1.0*CountDic.get('A',0)/Tot
                                  RecBacSup.append('\t'.join([RB,str(M),N,str(CountDic.get('T',0)),str(CountDic.get('A',0)),str(Mpro),str(Wpro)]))
                                  if RB=='B':
                                       if Wpro>FN:
                                           NodeID=Pos2NodeID[MostDec]

                                           MutIn=fiNodeMut[NodeID]
                                        
                                           MutIn.append('Back '+str(M))
                                           fiNodeMut[NodeID]=MutIn
                                                                                  
                                           fiNode2In[NodeID]=NodeID+' [label=\"'+NodeID+'\\n'+'\\n'.join(MutIn)+'\"]\n'
                                        #   print ('back added',M,NodeID,fiNode2In[NodeID])
                                  elif RB=='R':
                                    #   print (Mpro,FP, Wpro,FN)
                                       if Mpro>FP:
                                           NodeID=Pos2NodeID[MostDec]

                                           MutIn=fiNodeMut[NodeID]
                                          
                                           MutIn.append('Rec '+str(M))
                                           fiNodeMut[NodeID]=MutIn
                                                                                 
                                           fiNode2In[NodeID]=NodeID+' [label=\"'+NodeID+'\\n'+'\\n'.join(MutIn)+'\"]\n'
                                         #  print ('rec added',M,NodeID,fiNode2In[NodeID]) 
    TarNode2AddFor=InvertDic1(BM2For)
  #  print (TarNode2AddFor)
    for TarNode in TarNode2AddFor:
        MutLs=TarNode2AddFor[TarNode]  
        
        fiNodeMut[Pos2NodeID[TarNode]]+=map(str,MutLs)
                                         

    out=['digraph G {\n']
    for NodeID in fiNodeMut:# fiNode2In:
       
        out.append(NodeID+' [label=\"'+NodeID+'\\n'+'\\n'.join(fiNodeMut[NodeID])+'\"]\n')
    for E in fiEdge2In:
        out.append(E+';\n')
    out.append('}')
    GetOut(OutGV,''.join(out))  
    GetOut(OutGV.replace(OutGV.split(os.sep)[-1],'')+'RecBackSupp.txt','\n'.join(RecBacSup))  
def TestBackRec2(BackRecGV,FinalGV,CellAnno,OriFas,FP,FN):
    OutGV=FinalGV[:-3]+'_backrec.gv'
    rbDec2Anc,rbNodeMut,rbNode2In,rbEdge2In	=ReadGV1(BackRecGV)
    rbMut2Node=InvertDic3(rbNodeMut)
#    print (rbMut2Node)
  
    fiDec2Anc,fiNodeMut0,fiNode2In0,fiEdge2In	=ReadGV1(FinalGV)
    fiNodeMut={}
    for i in fiNodeMut0:
        In=fiNodeMut0[i]
        New=[]
        for ii in In:
            if ii[0]!='R' and ii[0]!='B' and ii[-1]!='A' and ii.find('?')==-1: New.append(ii)
        fiNodeMut[i]=New
      
    fiNode2In={}
    for NodeID in fiNode2In0:
        MutIn=fiNodeMut[NodeID] 
              
        fiNode2In[NodeID]=NodeID+' [label=\"'+NodeID+'\\n'+'\\n'.join(MutIn)+'\"]\n'        
     
    Pos2NodeID={}
    for fD in fiDec2Anc:
        ID=int(fD.split('Node')[-1])
        Pos2NodeID[ID]=fD
    NodeID2Pos=InvertDic(Pos2NodeID)    
    
    TH2CellLs=ReadTHAnno1(CellAnno,1,0)
   
    if OriFas[-3:]=='meg': CellLs,Cell2Seq=ReadMegSeq(OriFas)
    else: CellLs,Cell2Seq=ReadFasSeq(OriFas)
    RecBacSup=['\t'.join(['Type','Mutation postion (from 0)','NodeID','Mutant Cell Count','Wild Cell Count','Mutant cell proportion','Mutant cell proportion'])]
    BM2For={}
    for N in rbNodeMut:
       MutLs=rbNodeMut[N]
       RecBaDic=GetRekBack(MutLs)
      
       for RB in RecBaDic:
          if RB!='F':
           MposLs=RecBaDic[RB]
           if len(MposLs)>0:
             for M in MposLs:
               ConFLs=RecBaDic['F']
               if len(ConFLs)>0:
                   AllDec=[]
               
                   LenLs=[]
                   Len2F={}
                   F2DecLs={}
                   for F in ConFLs:
                     if F in Pos2NodeID:
                      Fid=Pos2NodeID[F]
                      DecLs=GetDecCla(Fid,fiDec2Anc)

                      CountTar=CountTarget(DecLs,fiNodeMut0,M)
                      print (RB,DecLs,MposLs,CountTar,ConFLs)
                      if CountTar>1: open('a','r').readlines()                      
                      LenLs.append(len(DecLs))
                      Len2F[len(DecLs)]=Len2F.get(len(DecLs),[])+[F]
                      F2DecLs[F]=DecLs
                   
                      AllDec+=DecLs
                   LenLs.sort()
               
                   AllDec=list(set(AllDec))
                   if len(LenLs)>0:
                     MostAnc=Len2F[LenLs[-1]][0]
                     MostDec=Len2F[LenLs[0]][0]
                
                     if len(AllDec)<LenLs[-1]: pass#print ('mutation order changed, so recurrent/back mutation cannot be added',MposLs)
                     else: 
                       DecCellLs=[]
                       for TH in F2DecLs[MostAnc]:
                           if TH in TH2CellLs: DecCellLs+=TH2CellLs[TH]
                   
                       Go='y'
                       if RB=='B': 
                            
                               Ancb=fiDec2Anc.get(Pos2NodeID[MostAnc],'')
                               AllAnc=[] #from dec
                               if Ancb!='':
                                   if Ancb!='root': AllAnc.append(NodeID2Pos[Ancb])
                                   while Ancb in fiDec2Anc:
                                         Ancb=fiDec2Anc[Ancb]
                                         if Ancb!='root': AllAnc.append(NodeID2Pos[Ancb])
                              
                               FirstConcFor=GetFirstConcMut(rbNodeMut[rbMut2Node['A'+str(M)+'T'][0]],AllAnc)    
                               if FirstConcFor!='': BM2For['A'+str(M)+'T']=FirstConcFor 
                               else: Go='n'
                                                         
                       if Go=='y':
                               CountDic=CountMut(DecCellLs,Cell2Seq,M) 
                               Tot=CountDic.get('T',0)+CountDic.get('A',0)
                              # print ('\n',N,RB,M,CountDic)                                
                               if Tot>0:
                                  Mpro=1.0*CountDic.get('T',0)/Tot
                                  Wpro=1.0*CountDic.get('A',0)/Tot
                                  RecBacSup.append('\t'.join([RB,str(M),N,str(CountDic.get('T',0)),str(CountDic.get('A',0)),str(Mpro),str(Wpro)]))
                                  if RB=='B':
                                       if Wpro>FN:
                                           NodeID=Pos2NodeID[MostAnc]

                                           MutIn=fiNodeMut[NodeID]
                                     
                                           MutIn.append('Back '+str(M))
                                           fiNodeMut[NodeID]=MutIn
                                                                             
                                           fiNode2In[NodeID]=NodeID+' [label=\"'+NodeID+'\\n'+'\\n'.join(MutIn)+'\"]\n'
                                          # print ('back added',M,NodeID,fiNode2In[NodeID])
                                  elif RB=='R':
                                     #  print (Mpro,FP, Wpro,FN)
                                       if Mpro>FP:
                                           NodeID=Pos2NodeID[MostAnc]

                                           MutIn=fiNodeMut[NodeID]
                                       
                                           MutIn.append('Rec '+str(M))
                                                                                
                                           fiNode2In[NodeID]=NodeID+' [label=\"'+NodeID+'\\n'+'\\n'.join(MutIn)+'\"]\n'
                                        #   print ('rec added',M,NodeID,fiNode2In[NodeID]) 
    TarNode2AddFor=InvertDic1(BM2For)
  # print (TarNode2AddFor)
    for TarNode in TarNode2AddFor:
        MutLs=TarNode2AddFor[TarNode]  
      
        fiNodeMut[Pos2NodeID[TarNode]]+=map(str,MutLs)
                                         

    out=['digraph G {\n']
    for NodeID in fiNodeMut:# fiNode2In:
      
        out.append(NodeID+' [label=\"'+NodeID+'\\n'+'\\n'.join(fiNodeMut[NodeID])+'\"]\n')
    for E in fiEdge2In:
        out.append(E+';\n')
    out.append('}')
    GetOut(OutGV,''.join(out))  
    GetOut(OutGV.replace(OutGV.split(os.sep)[-1],'')+'RecBackSupp.txt','\n'.join(RecBacSup))  
def CountTarget(DecLs,NodeMut,T):
    C=0
    for D in DecLs:
        MutLs=NodeMut[D]
        for M in MutLs:
            M=M.replace('Rec','').replace('Back','').replace(' ','').replace('A','').replace('T','')
            if int(M)==T:  C+=1
    return C            
def GetFirstConcMut(ConcMutLs,AllAnc):
   
    AllAnc.reverse()    
   
    for i in AllAnc:
        if 'A'+str(i)+'T' in ConcMutLs: 
           First=i
           return i
  #  print ('no concurrent mutation')
    return ''    
           
def CountMut(SeqLs,SeqDic,Pos):
    Dic={}
    for Seq in SeqLs:
        if '>'+Seq in SeqDic: Nuc=SeqDic['>'+Seq][Pos]
        else: Nuc=SeqDic['#'+Seq][Pos]
        Dic[Nuc]=Dic.get(Nuc,0)+1
    return Dic        
                   
    
def GetRekBack(MutLs):
    Dic={'R':[],'B':[],'F':[]}
    for M in MutLs:
      if M.find('?')==-1:  
        if M.find('Rec')!=-1: Dic['R'].append(int(M.replace('Rec','').replace('A','').replace('T','')))
        elif M[-1]=='A':  Dic['B'].append(int(M.replace('Rec','').replace('A','').replace('T','')))  
        else: Dic['F'].append(int(M.replace('Rec','').replace('A','').replace('T','')))          
    return Dic
def SupportRecBack(Clo2Cell,Cell2Seq,Type2Pos,OutF):
    out=['\t'.join(['Type','Position from 0','Clone','Cell count Wild','Cell count Mut','Cell count Miss'])+'\n']
    for Type in Type2Pos:
        PosLs=list(set(Type2Pos[Type]))
        for Pos in PosLs:
            for Clo in Clo2Cell:
                CellLs=Clo2Cell[Clo]
                Dic={}
                for Cell in CellLs:
                    if '>'+Cell in Cell2Seq:
                         Nuc=Cell2Seq['>'+Cell][Pos]
                    else:  Nuc=Cell2Seq['#'+Cell][Pos]
                    Dic[Nuc]=Dic.get(Nuc,0)+1
                out.append('\t'.join([Type,str(Pos),Clo,str(Dic.get('A',0)),str(Dic.get('T',0)),str(Dic.get('?',0))])+'\n')
    GetOut(OutF,''.join(out))                

def GetAncTilFor(Bmut,Bnode,Dec2Anc,NodeMut):
    ForMut='A'+Bmut[1:][:-1]+'T'
  #  print (Bnode,Bmut,ForMut)
    AncLs=[Bnode]	
    Find='n'	
    while Find=='n':
      if Bnode not in Dec2Anc: Find='y'	
      else:	  
       Bnode=Dec2Anc[Bnode]
       AncLs.append(Bnode)	   
       Mls=NodeMut[Bnode]	
        
       if Mls.count(ForMut)!=0: Find='y'
    ForN=Bnode	   
  
    AllDec=Getalldec(ForN,Dec2Anc)
  
    ONls=[]
    for A in AllDec:
       if AncLs.count(A)==0: ONls.append(A)
  	
    return ForN,AncLs,ONls,GetMut(AncLs,NodeMut),GetMut(ONls,NodeMut)	
def GetMut(Nls,NodeMut):
    Mls=[]
    for N in Nls:
        Mls+=NodeMut[N]
    return Mls	
def MapCOI(Chi,Patls,Chi2Pat2COI):
    Cls=[]
    for P in Patls:
       if P.strip()!='':
          if Chi[1:] in Chi2Pat2COI:
              Ci0=Chi2Pat2COI[Chi[1:]]  
              if P[1:] in Ci0:              
                  Cls.append(float(Chi2Pat2COI[Chi[1:]][P[1:]]))
    return Cls		
def AnaBack(GV,COI):
   
    Out=GV[:-3]+'_back.gv'
    
    Dec2Anc,NodeMut,Node2In,Edge2In	=ReadGV(GV)
    Dec2Anc2COI=ReadCOI(COI)
    
    BackM2NodeLs={}
    for N in NodeMut:
        MutLs=NodeMut[N]
        for M in MutLs:
          if M.strip()!='':	
           if M[-1]=='A': BackM2NodeLs[M]=BackM2NodeLs.get(M,[])+[N]
 #   print ('back list',BackM2NodeLs)
    	 
    	  
    Node2RmBack={}
    for B in BackM2NodeLs:
        NodeLs=BackM2NodeLs[B]
        for N in NodeLs:
            
            ForN,AncNLs,OthNLs,ForToBackMls,ForToNonBackMls=GetAncTilFor(B,N,Dec2Anc,NodeMut)			
            ForToBackMls.remove(B)
            if 'A'+B[1:][:-1]+'T' in ForToBackMls: ForToBackMls.remove('A'+B[1:][:-1]+'T')		
           # print (B,N,ForToBackMls,ForToNonBackMls)
            SuppCOI=MapCOI(B,ForToBackMls,Dec2Anc2COI)
            AltCOI=MapCOI(B,ForToNonBackMls,Dec2Anc2COI)
            P=999		
            if len(SuppCOI)>0 and len(AltCOI)>0:
                aveS=sum(SuppCOI)/len(SuppCOI)
                aveA=sum(AltCOI)/len(AltCOI)
               		
                if aveA<aveS:
                   A=stats.ttest_ind(SuppCOI,AltCOI)			
                   P=A.pvalue
            elif len(SuppCOI)>0:
                 aveS=sum(SuppCOI)/len(SuppCOI)
                 if aveS>0.1:			 
                     A=stats.ttest_1samp(SuppCOI,0.1)
                     P=A.pvalue	
            elif len(AltCOI)>0:
                 aveA=sum(AltCOI)/len(AltCOI)
                 if aveA<0.1:			 
                     A=stats.ttest_1samp(AltCOI,0.1)
                     P=A.pvalue				 
          
            if P>0.05:
               Node2RmBack[N]=Node2RmBack.get(N,[])+[B]	
 #   print ('remove bad backward',Node2RmBack)		   
    out='digraph D {\n'
    for N in Node2In:
        if N not in Node2RmBack: out+=Node2In[N]
        else:
           In=Node2In[N]
           MutLsO=NodeMut[N]
           RmLs=Node2RmBack[N]
           NewIn=[]
           for OM in MutLsO:
                if RmLs.count(OM)==0: NewIn.append(OM)
                else:
                    print ('back mutation removed',N,OM)								
           out+=N+'  [label=\"'+N+'\\nMut:'+';'.join(NewIn)+'\\n\"]\n'	   
    for E in Edge2In:
        out+=E+'\n'
    out+='}'
    GetOut(Out,out)	
def CompMutFre(OriFas): #OriFas_MutFre.txt
    if OriFas[-4:]=='.meg': 
        CellLs,Cell2Seq=ReadMegSeq(OriFas)
        Out=OriFas[:-4]+'_MutFre.txt' 
    else: 
        CellLs,Cell2Seq=ReadFasSeq(OriFas)
        Out=OriFas[:-6]+'_MutFre.txt'  
    Pos2C={}
    Len=len(Cell2Seq[CellLs[0]])
    for Cell in CellLs:
        Seq=Cell2Seq[Cell]
        if Seq.find('T')!=-1:
            c=0
            while c<Len:
                Nuc=Seq[c]
                if c not in Pos2C: Pos2C[c]={}
                Pos2C[c][Nuc]=Pos2C[c].get(Nuc,0)+1
                c+=1
    Pos2Fre={}
    out='Pos from 0\tMut\tWild\tMiss\tMutFre\n'
    c=0
    while c<Len:
        Cou=Pos2C[c]
        M=Cou.get('T',0)
        W=Cou.get('A',0)
        Miss=Cou.get('?',0)
        Tot=M+W
        if Tot==0: Fre=0
        else:
            Fre=1.0*M/Tot 
        Pos2Fre[c]=Fre
        out+=str(c)+'\t'+'\t'.join(map(str,[M,W,Miss,Fre]))+'\n'
        c+=1
    GetOut(Out,out)
   
    return Pos2Fre   
def GetPosMut(MutIDLab,MutIDLs): 
    MutID=MutIDLab.replace('Rec','').replace('Back','').replace('A','').replace('T','').replace('?','').strip()
    if MutIDLs.count(MutID)!=0: Pos=MutIDLs.index(MutID)
    else: Pos= MutIDLs.index(MutID+'\n')
   
    return Pos    
def TopHapGV2MOAinGV(GVin,Pos2MutFre,MutIDLs): #_MutTree.txt
     Dec2Anc,NodeMut,Node2In,Edge2In=ReadGV1(GVin)
    
     AllMutLs=[]
     BackRecLs=[]
     EdgeIn=[]
     NodeIn=[]
     MissEdge={}
     Node2FirstMut={}
     
     Node2LastMut={}
     MutTreeIn=['\t'.join(['Branch','Ancestor','Muttation','Loss of mutation','Parallel mutations','Cell count'])]
     for E in Edge2In:
        if E.find('->')!=-1:
            AncDec=E.replace(';','').strip().split('->')

            for N in AncDec:
                 MutIn=NodeMut[N]
                 for M in MutIn:
                     Pos=GetPosMut(M,MutIDLs)
                     if M[0]=='A' and M[-1]=='T': pass
                     else:  
                         BackRecLs.append(Pos)
     for E in Edge2In:
        if E.find('->')!=-1:
            AncDec=E.replace(';','').strip().split('->')
            Last=''
            Anc='y'
            First=''
            for N in AncDec:
                 MutIn=NodeMut[N]
                 Mut2Freq={}
                 OtherMut={'B':[],'R':[],'Amb':[]}
                 for M in MutIn:
                     Pos=GetPosMut(M,MutIDLs)
                     if M[0]=='A' and M[-1]=='T':
                         AllMutLs.append(Pos)
                         Mut2Freq[Pos]=Pos2MutFre[Pos]
                     else:  
                      
                         if M.find('?') !=-1:  OtherMut['Amb'].append(Pos)                        
                         elif M[-1]=='A' :  OtherMut['B'].append(Pos) 
                         elif M[-1]=='T' :  OtherMut['R'].append(Pos) 
                         else: open('a','r').readlines()                         
                 Freq2Mut=InvertDic1(Mut2Freq)  
                 FreqLs=list(Freq2Mut.keys())
                 FreqLs.sort(reverse=True)
               
                 MOr=[]
                 for Fre in FreqLs:
                     MutLs0=Freq2Mut[Fre]
                     MutLs=[]
                     for m0 in MutLs0:
                         if m0 not in BackRecLs: 
                             MutLs.append(m0)
                                                  
                     MOr+=list(map(str,MutLs))
                
                 if Anc=='y' and len(MOr)>0 and First=='': First=MOr[0]
                 if Anc=='y' and N not in Dec2Anc:
                     MutTreeIn.append('\t'.join([N,'root',';'.join(MOr),';'.join(map(str,OtherMut['B'])),';'.join(map(str,OtherMut['R'])),';'.join(map(str,OtherMut['Amb']))]))
                 elif Anc!='y':
                     MutTreeIn.append('\t'.join([N,AncDec[0],';'.join(MOr),';'.join(map(str,OtherMut['B'])),';'.join(map(str,OtherMut['R'])),';'.join(map(str,OtherMut['Amb']))]))                 
                 if Anc=='y' and len(MOr)!=0: Last=MOr[-1]   
                               
                 Len0=len(MOr)-1
                 if len(MOr)!=0:Node2FirstMut[N]= MOr[0]
                 c=0
                 while c<Len0:
                        EdgeIn.append('\"'+MOr[c]+'\"->\"'+MOr[c+1]+'\";\n')
                        NodeIn+=['\"'+MOr[c]+'\" [label=\"'+MOr[c]+'\\nMut:'+MOr[c]+'\"];\n','\"'+MOr[c+1]+'\" [label=\"'+MOr[c+1]+'\\nMut:'+MOr[c+1]+'\"];\n']                        
                        c+=1

                 Anc='n'
                 if len(MOr)!=0: Node2LastMut[N]=MOr[-1]
        if Last!='' and len(MOr)!=0:
                        EdgeIn.append('\"'+Last+'\"->\"'+MOr[0]+'\";\n')
                     
        if First=='' and AncDec[0]!='root':                 
         
           if len(MOr)>0: MissEdge[AncDec[0]]=MissEdge.get(AncDec[0],[])+[MOr[0]]
       # print (AncDec,MOr,First)
            
  #   print ('mising edge',MissEdge)
   
     for Miss in MissEdge:
        Find='n'
        Anc=Miss
        while Find=='n':
           Anc=Dec2Anc[Anc]
           Find=Node2LastMut.get(Anc,'n')
           if Anc=='root':
               Find=Node2LastMut.get(Anc,'root')
        Ls=MissEdge[Miss]
        for i in Ls:        
           EdgeIn.append('\"'+Find+'\"->\"'+i+'\";\n')   
           print ('add \"'+Find+'\"->\"'+i+'\";\n')
                
     GetOut(GVin[:-3]+'_MutTree.txt','\n'.join(MutTreeIn))
     
     Root=''
     A=list(Dec2Anc.keys())[0]

        
     while A in Dec2Anc:
         if Dec2Anc[A]=='root': break
         A=Dec2Anc[A]
         
     if 'root' in Node2FirstMut:  Root=Node2FirstMut['root']    
     else: 
        if A in Node2FirstMut: Root=Node2FirstMut[A] 
  
     if Root!='': EdgeIn.append('root->\"'+Root+'\";\n')     
     out='digraph D {\n'+''.join(list(set(NodeIn)))+''.join(list(set(EdgeIn)))+'}\n'
     GetOut(GVin[:-3]+'_MOA.gv',out) 
  
     Miss2Fre={}
     for P in Pos2MutFre:
        if P not in AllMutLs and P not in BackRecLs: 
            Miss2Fre[P]= Pos2MutFre[P]       
    # print (len(Miss2Fre))  
     GoodLs=[]
     AllMutLs=list(set(AllMutLs))
     for i in AllMutLs:
         if i not in BackRecLs: GoodLs.append(i)     
     return Miss2Fre,GoodLs
def MOAoutFixLabel(MOAoutgv): #1.txt'
    Dec2Anc,NodeMut,Node2In,Edge2In	=ReadGVmoa(MOAoutgv)   
    out=['digraph G {\n']
    for Node in NodeMut:
        Label=NodeMut[Node][0]
        out.append('\"'+Label+'\"'+' [label=\"'+Label+'_'+Node+'\"];\n')
    for E in Edge2In:
        E=E.split('->')
        out.append('\"'+NodeMut[E[0].strip()][0]+'\"'+' -> '+'\"'+NodeMut[E[1].strip()][0]+'\"'+';\n')
    out.append('}\n')
    GetOut(MOAoutgv[:-4]+'1.txt',''.join(out))    
def MakeMOAin2(PosLs,Miss,OriFas):
    if OriFas[-3:]=='meg': 
        CellLs,Cell2Seq=ReadMegSeq(OriFas)
        OutFa=OriFas[:-4]+'MOAIn.fasta'
        Out=OriFas[:-4]+'MOAlabelIn.txt'
    else:
        CellLs,Cell2Seq=ReadFasSeq(OriFas)
        OutFa=OriFas[:-6]+'MOAIn.fasta'
        Out=OriFas[:-6]+'MOAlabelIn.txt'        
    out='label\tdays\n'
    if Miss!='': PosIDLs=PosLs+[Miss]
    else: PosIDLs=PosLs

    for i in PosIDLs:
        out+=str(i)+'\t0\n'

    GetOut(Out,out)    
    out=[]
    for Cell in CellLs:
        Seq=Cell2Seq[Cell]
        Seq1=''
      
        for P in PosIDLs:
               
                Seq1+=Seq[P]
        if Seq1.find('T')!=-1:        
            out+=[Cell.replace('#','>')+'\n'+Seq1+'\n']    
    GetOut(OutFa,''.join(out))
    return len(out) ,Out,OutFa   
def MakeMOAin1(PosID,GV,Fas):
 
    Dec2Anc,NodeMut,Node2In,Edge2In	=ReadGV(GV)
    CellLs,Cell2Seq=ReadFasSeq(Fas)
    Out=PosID[:-4]+'MOAlabelIn.txt'
    OutFa=Fas[:-6]+'MOAIn.fasta'
    out='label\tdays\n'
    PosID=open(PosID,'r').readlines()
    PosIDLs=[]
    for i in PosID:
        out+=i.strip()+'T\t0\n'
        PosIDLs.append(int(i))	
    BackLs=[]
    for N in NodeMut:
        Mls=NodeMut[N]
        for M in Mls:
            if M[0]=='T': BackLs.append(int(M[1:][:-1]))
    BackLs=list(set(BackLs))
    print (BackLs)
    NucDic={'A':'T','T':'A','?':'?'}
    for B in BackLs:
        out+=str(B)+'A\t0\n'
        for Cell in CellLs:
            Nuc=Cell2Seq[Cell][PosIDLs.index(B)]
            NucO=NucDic[Nuc]		
        
            Cell2Seq[Cell]+=NucO		
 
    GetOut(Out,out)    
    out=''
    for Cell in CellLs:
        out+=Cell+'\n'+Cell2Seq[Cell]+'\n'
    GetOut(OutFa,out)
def FillChaMat(File,Min):
    Site2In=ListCol(File)
    Site2Keep={}
    for Site in Site2In:
       if Site!='cellBC':
        In=Site2In[Site]
        Cha2C={}
        for i in In:
           Cha2C[i]=Cha2C.get(i,0)+1
        for Cha in Cha2C:
            if Cha2C[Cha]>=Min: Site2Keep[Site]=Site2Keep.get(Site,[])+[Cha]
            else: print ('remove due to small freq',Site,Cha,Cha2C[Cha])
   # print (Site2Keep)
    SiteOr=list(Site2Keep.keys())
    out=['cellBC\t'+'\t'.join(SiteOr)+'\n']
    Len=len(Site2In['cellBC'])
    c=0
    while c<Len:
       In=[Site2In['cellBC'][c]]
       for Tar in SiteOr:
           Ind=Site2In[Tar][c]
           if Ind in Site2Keep[Tar] or Ind=='-' or Ind=='?': In.append(Ind)
           else: In.append('0')
       c+=1
       out.append('\t'.join(In)+'\n')
    OutF=File[:-4]+'_Fill'+str(Min)+'.txt'
    GetOut(OutF,''.join(out))
    return OutF    
        
def ChaMat2Fas(File):
   
    RmNote=[]
    if os.path.exists(File)!=True:
        print (File,'does not exist')
    else:	
        print (File)
        PosOut=File[:-4]+'_Position.txt'	
        Cha2Val=ListCol(File)	
        ID=File.split('\\')[-1].replace('_character_matrix.alleleThresh.txt','')
    
        File=open(File,'r').readlines()
        ChaOrder0=File[0].strip().split('\t')[1:]
        CellCol=File[0].strip().split('\t')[0]	
        Posout='Character\tState\tposition in Fas\n'
        ChaDic={}
        Pos=0
        ChaOrder=[]	
        for Cha in ChaOrder0:
            ValLs=list(set(Cha2Val[Cha]))
            ValLs.sort()
            MissC=Cha2Val[Cha].count('-')
            MissPeo=1.0*MissC/len(Cha2Val[Cha])
         	
            if MissPeo>=0.3: 
               # print ('remove',Cha,'due to >=30% missing data',MissC,len(ValLs),MissPeo)
                RmNote.append(Cha+'\t'+str(MissC)+'\t'+str(len(Cha2Val[Cha]))+'\t'+str(MissPeo)+'\n')			
            else:
             ChaOrder.append(Cha)		
    			
             ChaDic[Cha]={}	
             StC=0		
             for Val in ValLs:
              if Val!='0' and Val!='-':		
                ChaDic[Cha][Val]=Pos	
                Posout+=Cha+'\t'+Val+'\t'+str(Pos+1)+'\n'
                Pos+=1
                StC+=1			
    		
        GetOut(PosOut,Posout)			
     #   print (ChaOrder)	
    
        File=File[1:]	
        Len=len(Cha2Val[ChaOrder[0]])
        C=0
        Fasout='>Normal\n'+('A'*Pos)+'\n'	
        while C<Len:	
           NewID=Cha2Val[CellCol][C]
           if str(C).find('000')!=-1 and str(C)[-1]=='0': print (NewID,C,Len)	   
    
           Seq=''	
           	   
           for Cha in ChaOrder:
               Cha2Pos=ChaDic[Cha]	   
               State=Cha2Val[Cha][C]
               PosLs=list(Cha2Pos)
               PosLs.sort()
    	   
               if State=='-': Seq+='?'*len(Cha2Pos)	
               elif State=='0': Seq+='A'*len(Cha2Pos)
               else:
                  Find='n'		   
                  for Pos in PosLs:
                      if Pos==State: 
                          Seq+='T'
                          Find='y'					  
                      else: Seq+='A'				  
                  if Find=='n':
    
                      open('A','r').readlines()	
           MissC=Seq.count('?')
           MissPro=1.0*MissC/len(Seq)
           if MissPro>=0.3:
    
               RmNote.append(NewID+'\t'+str(MissC)+'\t'+str(len(Seq))+'\t'+str(MissPro)+'\n')			   
           else:		   
               Fasout+='>'+NewID+'\n'+Seq+'\n'				  
           C+=1	
        	   
        GetOut(PosOut[:-4]+'.fasta',Fasout)
        GetOut(PosOut[:-4]+'_Rm.txt','Site/Cell\tMissing count\tTotal count\tProportion\n'+''.join(RmNote))	
def Fas2Json2(Fa):
 
    CellLs,Cell2Seq=ReadFasSeq(Fa) 
    print ('outgroup sequence assigned: ',CellLs[0])
    OutSeqID=CellLs[0]
    OutSeq=Cell2Seq[OutSeqID]
    Len=len(OutSeq)
    GetOut(Fa[:-3]+'OutSeq.fasta',OutSeq)	
    In=[]
    for Cell in CellLs:
          Seq=Cell2Seq[Cell]
          if Seq==OutSeq:		
            print ('skiped because it is outgroup sequence: ',Cell,len(Seq))
          else:		
            VarLs=[]
            C=0
            while C<Len:
                if Seq[C]!=OutSeq[C]: VarLs+=[str(C),'\"'+Seq[C]+'\"']
                C+=1
            In.append('{'+'\"V\":['+','.join(VarLs)+'],\"I\":[\"20220101\",\"'+Cell.replace('>','')+'\",\"NA\",\"'+Cell.replace('>','')+'\",\"'+Cell.replace('>','')+'\"]}')
    out='{\"sequences\":['+','.join(In)+']}'
    GetOut(Fa[:-3]+'.json',out)	#count from 0
    	   
    
def MakeAncSeq(Ta,Ref):

    NodeMap=Ta.replace('.csv','_nodeMap.txt')

    Out=Ta.replace('.csv','WithAnc.meg')
    NodeMap=open(NodeMap,'r').readlines()
    NodeMap=NodeMap[1:]
    Node2Label={}
    RefNodeID=''
    Dec2Anc={}
    Label2Count={}
    Label2Seq={}
    Label2Cell={}
    
    for i in NodeMap:
        i=i.strip().split('\t')
        ii=[]
        for Item in i:
            Item=Item.strip()
            if Item!='': ii.append(Item)
        Label=ii[0].strip().replace(' ','_')
        Node=ii[1].strip()
        Decs=[ii[2].strip(),ii[3].strip()]
        if Label[:5]=='Node_': 
            for Dec in Decs:
                Dec2Anc[Dec]=Node
        Label2Cell['Node_'+Node]=Label			
        Label2Count[Label]=0
        Label2Seq['Node_'+Node]=''
        if Label==Ref or Label=='Normal': RefNodeID=Node
        Node2Label[Node]=Label

    RefAnc=Dec2Anc[RefNodeID]
    Code=RefAnc in Dec2Anc
    RmNode=Node2Label[RefAnc]
    if Code==True: print ('Error: Root is incorrect. Fix the nwk file and redo ancestor analysis', Ta)
    else:
        AddNode=[]
        Ta=open(Ta,'r').readlines()
        Head=Ta[0].strip().split(',')
        Node2Col={}
        Node2Prob={}
        c=0
        Len=len(Head)
        while c<Len:
            i=Head[c].strip()
            Code=i in Label2Seq
            if Code==True: Node2Col[i]=c
            c+=1
        Ta=Ta[1:]
    
        PosiP=1
        for i in Ta:
            i=i.strip().split(',')
            PosiC=int(i[0])
            Nuc=i[1]
            if PosiC!=PosiP:
                PosiP=PosiC
                Node2BestNuc=GetBestNuc(Node2Prob)
                for Node in Node2BestNuc:
                    Label2Seq[Node]+=Node2BestNuc[Node]
                Node2Prob={}
            for Node in Node2Col:
                Prob=i[Node2Col[Node]]
                Code=Node in Node2Prob
                if Code!=True: Node2Prob[Node]={}
                Node2Prob[Node][Nuc]=float(Prob)
    
    Node2BestNuc=GetBestNuc(Node2Prob)
    for Node in Node2BestNuc:
                    Label2Seq[Node]+=Node2BestNuc[Node]

    out=''
  
    for Label in Label2Seq:

     if Label2Seq[Label]!='':
          out+='#'+Label+'CellID'+Label2Cell[Label]+'\n'+Label2Seq[Label]+'\n'

    GetOut(Out,out)      
def GetOut(OutFile,outIn):
 OutF=open(OutFile,'w')
 OutF.write(outIn)
 OutF.close()
def GetFasSeq(Dic, OutF):
   out=[]
   for i in Dic:
       S=Dic[i]
       if i[0]!='>': i='>'+i.replace('#','')
       out.append(i+'\n'+S+'\n')
   GetOut(OutF,''.join(out))    