import os
import shutil
import sys
import Functions3

HapOr=sys.argv[3]
OriFasHap=sys.argv[2] #for MOA
OriFas=sys.argv[4]
COIcut=float(sys.argv[5])
COImatrix=sys.argv[1]#'C:\\Users\\kumarlab\\Desktop\\FastMOA\\SingleCell\\COIfromOriFas\\TallG5clone_40_1_BEAMin_BEAM.faHap\\TopHapOut\\MOAout\\COI_matrix.txt'
GVin=COImatrix.replace('MOAout'+os.sep+'COI_matrix.txt','TopHap_MutTree.gv')
prune='s'
c=0

print (GVin)
Functions3.MapCOI3(COImatrix,GVin)
print ('prune unsupported mut by COI')
Functions3.SumCOI(COImatrix,GVin[:-3]+'_COI3.gv',COIcut) #_COI3_ave3.gv'

def CopyAndDel(File,New):	
   shutil.copy(File, New)
   os.remove(File)	

if GVin.find('_prune.gv')==-1:
    shutil.copy(GVin,GVin[:-3]+'_prune.gv') 
    shutil.copy(GVin[:-3]+'_COI3.gv',GVin[:-3]+'_prune_COI3.gv')		
    shutil.copy(GVin[:-3]+'_COI3_ave3.gv',GVin[:-3]+'_prune_COI3_ave3.gv')	
elif GVin==GVin.replace(GVin.split('\\')[-1],'')+'TopHap_MutTree_prune.gv': pass
else:
   CopyAndDel(GVin,GVin.replace(GVin.split('\\')[-1],'')+'TopHap_MutTree_prune.gv')
   CopyAndDel(GVin[:-3]+'_COI3.gv',GVin.replace(GVin.split('\\')[-1],'')+'TopHap_MutTree_prune_COI3.gv') #TopHap_MutTree_prune_COI3.gv
   CopyAndDel(GVin[:-3]+'_COI3_ave3.gv',GVin.replace(GVin.split('\\')[-1],'')+'TopHap_MutTree_prune_COI3_ave3.gv') #TopHap_MutTree_prune_COI3_ave3.gv
   CopyAndDel(GVin[:-3]+'_COI3_ave3.png',GVin.replace(GVin.split('\\')[-1],'')+'TopHap_MutTree_prune_COI3_ave3.png') #TopHap_MutTree_prune_COI3_ave3.png
   CopyAndDel(GVin[:-3]+'_COI3_AveCOI3.txt',GVin.replace(GVin.split('\\')[-1],'')+'TopHap_MutTree_prune_COI3_AveCOI3.txt') #TopHap_MutTree_prune_COI3_AveCOI3.txt

print ('test recurrent and back')
FP,FN, TH2Count=Functions3.CompFPFN(OriFasHap[:-11]+'.fasta',GVin.replace(GVin.split('\\')[-1],'')+'TopHapNA__2022_1_FillAllCell_anno1.txt')
print ('FP,FN',FP,FN)

GVforTest=GVin.replace(GVin.split('\\')[-1],'')+'TopHap_MutTree_prune_COI3_ave3.gv'
Functions3.AnaRecBack(GVforTest,COImatrix,OriFas,FP,FN)
Functions3.GetOut(GVin.replace(GVin.split('\\')[-1],'')+'FPFN.txt','FP: '+str(FP)+'\n'+'FN: '+str(FN)+'\n')