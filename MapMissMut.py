import os
import sys
import Functions3

PosIDFile='NA'
OutFol=sys.argv[1]#'E:\\Desktop\\fastMOA\\SingleCell\\TopHapAll\\TallG6cell_500_5_BEAMin_BEAM.faHap\\TopHapOut\\'
GVin=sys.argv[2]#'E:\\Desktop\\fastMOA\\SingleCell\\TopHapAll\\TallG6cell_500_5_BEAMin_BEAM.faHap\\TopHapOut\\MutTree.gv'
OriFas=sys.argv[3]#'E:\\Desktop\\fastMOA\\SingleCell\\TopHapAll\\TallG6cell_500_5_BEAMin.meg'
BEAMOriFas=sys.argv[4]#'E:\\Desktop\\fastMOA\\SingleCell\\TopHapAll\\TallG6cell_500_5_BEAMin_BEAM.meg'
COICut=sys.argv[5]#'0.4'
MinCellC=int(sys.argv[6])#3 #>=3 include 
python=sys.argv[8]#'python'

MinCellC1=int(sys.argv[7])#2 #>=2 cells have lost mutation, then attach it
FPFNFile=OutFol+'FPFN.txt'
FPFN=open(FPFNFile,'r').readlines()
FP=float(FPFN[0].split(':')[1])#0.05049711080897349
FN=float(FPFN[1].split(':')[1])
print ('compute mutation freq')
Pos2MutFre=Functions3.CompMutFre(OriFas) #OriFas_MutFre.txt
SNVc=len(Pos2MutFre)
FreqTa=OriFas[:((-1)*(len(OriFas.split('.')[-1])+1))]+'_MutFre.txt'
print ('SNV count',SNVc)


if os.path.exists(PosIDFile)==True: PosIDOr=open(PosIDFile,'r').readlines()
else: 
    PosIDOr=list(map(str,range(0,len(Pos2MutFre)+1)))
 
MissMut2Freq,GoodPosLs=Functions3.TopHapGV2MOAinGV(GVin,Pos2MutFre,PosIDOr) #GVin_MOA.gv GVin_MutTree.txt

MutTreeFile=GVin[:-3]+'_MutTree.txt'

Freq2Mis=Functions3.InvertDic1(MissMut2Freq)  
FreqLs=list(Freq2Mis.keys())
FreqLs.sort(reverse=True)
#print (GoodPosLs,MissMut2Freq)

MOAupGV=GVin[:-3]+ '_MOA.gv'
c=0
Go='y'
for Freq in FreqLs:
    MisLs=Freq2Mis[Freq]
   # print (Freq,MisLs,GoodPosLs,len(GoodPosLs))
    for Mis in MisLs:
      if Go=='y':
       # print (Mis)
        SeqNum,LabFile,MOAin1=Functions3.MakeMOAin2(GoodPosLs,Mis,OriFas) #OriFas[:-4]+'MOAIn.fasta','MOAlabelIn.txt'
       # print (SeqNum)
        os.system(python+' FastMOA_update.py '+MOAin1+' '+LabFile+' --initial_tree '+ MOAupGV+' -o '+OutFol+'MOAmap'+str(Mis)+' --threshold '+COICut+' --disable_graph_flipping --flip_pass_thresh 1.1') #--lock_init_tree')
       # print (Mis, python+' FastMOA_update.py '+MOAin1+' '+LabFile+' --initial_tree '+ MOAupGV+' -o '+OutFol+'MOAmap'+str(Mis)+' --threshold '+COICut+' --disable_graph_flipping --flip_pass_thresh 1.1') #--lock_init_tree')
        MOAupGV0=OutFol+'MOAmap'+str(Mis)+os.sep+MOAin1.split(os.sep)[-1][:-6]+'_2.txt'
        if os.path.exists(MOAupGV0)==True:
            Functions3.MOAoutFixLabel(MOAupGV0) #1.txt'
            MOAupGV=OutFol+'MOAmap'+str(Mis)+os.sep+MOAin1.split(os.sep)[-1][:-6]+'_21.txt'
          #  print (MOAupGV)
            GoodPosLs.append(Mis)            
        else: Go='n'    


Dec2Anc,NodeMut,Node2In,Edge2In=Functions3.ReadGV1(MOAupGV)  

print ('filter single intermediate branching')

Dec2Anc,NodeMut,Node2In,Edge2In=Functions3.ReadGV1(MOAupGV)  

Anc2Dec=Functions3.InvertDic1(Dec2Anc)
ProMut=[]
for D in Dec2Anc:
    if D not in Anc2Dec:
        A=Dec2Anc[D]
        Dls=Anc2Dec[A]        
        if len(Dls)>1: ProMut.append(int(D.replace('\"','')))
print ('inter single attach',ProMut)   
AllMut=list(set(list(Dec2Anc.keys())+list(Anc2Dec.keys())))

GoodPosLs=[]
for M in AllMut:
    M=M.replace('\"','')
    if M!='root':
        if int(M) not in ProMut :  GoodPosLs.append(int(M))
#print (GoodPosLs,len(GoodPosLs))    

out=['digraph G {\n']
for N0 in Node2In:
    N=N0.replace('\"','')
    if N=='root': out.append(Node2In[N0])
    else:
        if int(N) not in ProMut: out.append(Node2In[N0])
for E0 in Edge2In:
 
    E=E0.replace('\"','').replace(';','').split('->')
    if E[0]=='root':
       if int(E[1]) not in ProMut: out.append(E0+';\n')
    else:
      if int(E[0]) not in ProMut and int(E[1]) not in ProMut: out.append(E0+';\n')
out.append('}\n')    

MOAupGV=OutFol+'MOAprune.gv'
Functions3.GetOut(MOAupGV,''.join(out))


MOAupGV=OutFol+'MOAprune.gv'
MOAupGV1=OutFol+'MOAprune.gv'

Go='y'
if Go=='y':
   print ('make expected and observed seq')
   MutP,ExtSeqDic=Functions3.GV2FasMoa(MOAupGV1,SNVc) 
   SeqNum,LabFile,MOAin1=Functions3.MakeMOAin2(MutP,'',OriFas) #OriFas[:-4]+'MOAIn.fasta','MOAlabelIn.txt'
   print ('annotate obs seq')
   if OriFas[-3:]=='meg':
       Functions3.SeqAnno(MOAupGV1[:-3]+'.fasta',OriFas[:-4]+'MOAIn.fasta','n') #MOAupGV[:-6]+'_CellAnnoAll.txt'
   else:
       Functions3.SeqAnno(MOAupGV1[:-3]+'.fasta',OriFas[:-6]+'MOAIn.fasta','n') #MOAupGV[:-6]+'_CellAnnoAll.txt'   
  
   AddedNodeLs=Functions3.AddBackTopHapMutTree(MOAupGV1[:-3]+'_CellAnnoAll.txt',MutTreeFile,MOAupGV1[:-3]+'.fasta',MOAupGV1,MinCellC,MinCellC1,FP,FN,SNVc,OriFas,MutP,'y') #MOAupGV_TopHap.gv fasta MOAupGV_TopHap0.gv
   

print ('test and remove')

ExpFas=MOAupGV1[:-3]+'_TopHap.fasta'
MOAupGV1=MOAupGV1[:-3]+'_TopHap0.gv'
if OriFas[-3:]=='meg': OriFasHap=OriFas[:-4]+'MOAIn.fasta'
else: OriFasHap=OriFas[:-6]+'MOAIn.fasta'

MutP,ExtSeqDic=Functions3.GV2FasMoa(MOAupGV1,SNVc) 
SeqNum,LabFile,MOAin1=Functions3.MakeMOAin2(MutP,'',OriFas)
Functions3.SeqAnno(ExpFas,OriFasHap,'n') #MOAupGV[:-6]+'_CellAnnoAll.txt' 
print ('test and remove intermediate')
BadIntermedLs=Functions3.TestIntermed(ExpFas,OriFas,FP,FN) #ExpFas[:-6] '_FilInter.fasta', need to remove inter clone from GV

UpExpFas=ExpFas[:-6] +'_FilInter.fasta'  

Functions3.SeqAnno(UpExpFas,OriFasHap,'n') #MOAupGV[:-6]+'_CellAnnoAll.txt' 
print ('test and remove doublet')
print ('\n\nIn Functions', UpExpFas[:-6]+'_CellAnnoAll.txt'  ,UpExpFas,OriFasHap,FP,FN)
Node2CellLs=Functions3.TestDoublet(UpExpFas[:-6]+'_CellAnnoAll.txt'  ,UpExpFas,OriFasHap,FP,FN) #final cell annotation

print ('final annotation')
UpExpFas=UpExpFas[:-6]+'_CellAnnoAll_final.fasta' 

Functions3.SeqAnno(UpExpFas,OriFasHap,'n') #MOAupGV[:-6]+'_CellAnnoAll.txt' 
print ('test and remove doublet 2')
Node2CellLs=Functions3.TestDoublet(UpExpFas[:-6]+'_CellAnnoAll.txt'  ,UpExpFas,OriFasHap,FP,FN) #final cell annotation

print ('position order: ',MOAupGV[:-3]+'_Pos.txt')
print ('MOA GV: ',MOAupGV)
print ('cell annotation: ',UpExpFas[:-6]+'CellAnnoAll_final_CellAnnoAll_final.txt')
print ('MOATopHapfasta: ',UpExpFas[:-6]+'CellAnnoAll_final_CellAnnoAll_final.fasta')
print ('MOATopHapGV: ',MOAupGV[:-3]+'_TopHap.gv')
print ('FP, FN',FP,FN)
print ('cell annotation before final: ',UpExpFas[:-6]+'CellAnnoAll_final.txt')
print ('MOATopHapfasta before final: ',UpExpFas[:-6]+'CellAnnoAll_final.fasta')