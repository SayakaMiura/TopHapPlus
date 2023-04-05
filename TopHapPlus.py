import sys
import os
import Functions3
import glob
import shutil
import glob


python='python'
COIcut=0.3
MOAcut=0.4
if sys.argv.count('-ST')!=0:
    COIcut=0.01
    MOAcut=0.1
    AllCellFas=sys.argv[sys.argv.index('-ST')+1]
    print (AllCellFas)   
if sys.argv.count('-COIcut')!=0: 
       COIcut=float(sys.argv[sys.argv.index('-COIcut')+1])
       print ('COI cutassigned',COIcut)
if sys.argv.count('-MOAcut')!=0: 
       MOAcut=float(sys.argv[sys.argv.index('-MOAcut')+1])
       print ('MOAcut assigned',MOAcut)   

print (COIcut,MOAcut)

MaxSNV=9999999999999 
MinVarError=5
def CopyAndDel(File,New):	
   shutil.copy(File, New)
   os.remove(File)	
def Clean(Tar):
   FileLs=glob.glob(Tar)
   for i in FileLs:
      os.remove(i)   
if sys.argv[1]=='-BEAM': 
    BEAMin=sys.argv[2]   
    Meg=BEAMin[:-4]+'_BEAM.meg'#sys.argv[2] #TallG5clone_10_1_BEAMin_BEAM.meg
    FinalRmLs=[BEAMin[:-4]+'_Hap.fasta',BEAMin[:-4]+'_HapMOAIn.fasta',BEAMin[:-4]+'_MutFre.txt',BEAMin[:-4]+'MOAIn.fasta',BEAMin[:-4]+'MOAlabelIn.txt',Meg[:-4]+'.fasta']
    if os.path.exists(Meg)!=True:
        CWD=os.getcwd()
        os.chdir('BEAM-master')
        os.system(python+' BEAM3.py '+BEAMin)
        os.chdir(CWD)
    VAF=sys.argv[3]	
    HF=sys.argv[4]
    BEAMinHap=BEAMin[:-4]+'_Hap.fasta'	
    DataID=BEAMin[:-4]	
    CellLs,Cell2Seq=Functions3.ReadMegSeq(Meg)
    if float(VAF)>1:
        VAF=str(1.0*float(VAF)/len(CellLs))
    if float(HF)>1:
        HF=str(1.0*float(HF)/len(CellLs))
    print (VAF,HF,len(CellLs))  
    In=[]
    for Cell in CellLs:
       In.append('>'+Cell.replace('#','')+'\n'+Cell2Seq[Cell]+'\n')	
    out='>Outgroup\n'+('A'*len(Cell2Seq[Cell]))+'\n'+''.join(In)
    Fas=Meg[:-4]+'.fasta'
    Functions3.GetOut(Fas,out)	
    Functions3.Fas2Json2(Fas)
    Json=Fas[:-3]+'.json'	
    if os.path.exists(Json)!=True: #count from 0
       print ('Json file was not created. Please check inout fasta file.')
    else:	
       OutSeq=Fas[:-3]+'OutSeq.fasta'	
       TopHapOut=Json[:-5]+'OutSeqTopHap.fasta'  
       Functions3.GetOut(TopHapOut,'>Outgroup\n'+open(OutSeq,'r').readlines()[0])   	   
       Hap=Json[:-5]+'Hap'	
       if os.path.exists(Hap)==True:
           print ('please delete ',Hap,' and try again')
       else:		   
           os.system(python+' vcf_json_parse.py '+Json+' --reference '+OutSeq+' --min_subgroup_size 3 --skip_mismatches --min_freq '+VAF+' -o '+Hap)
           Functions3.GetVAFpos(Hap+os.sep+'variant_stats.txt',float(VAF))	#Hap+os.sep+'variant_stats.txt'
           os.remove(Hap+os.sep+'Haplotypes.txt')
           os.remove(Hap+os.sep+'NA__2022_1.fasta')	
           os.remove(Hap+os.sep+'sequence_annotations.txt')
           os.remove(Hap+os.sep+'variant_stats.txt')		   
           os.system(python+' vcf_json_parse.py '+Json+' --reference '+OutSeq+' --min_subgroup_size 3 --skip_mismatches --min_freq 1.0 -o '+Hap+' --pos_file '+Hap+os.sep+'variant_stats_Force.txt') #position from 1			   
           print ('make haplotype for the original')
           CellInLs,CellIn2Seq=Functions3.ReadMegSeq(BEAMin)
           VarSite=open(Hap+os.sep+'Haplotypes.txt','r').readlines()
           out=''
           for C in CellInLs:
                Seq=CellIn2Seq[C]
                New=''				
                for V in VarSite:
                    New+=Seq[int(V)]
                out+='>'+C.replace('#','')+'\n'+New+'\n'         			
           Functions3.GetOut(BEAMinHap,out)				           		   
elif sys.argv[1]=='-Cas': 
    ChaMat0=sys.argv[2]
    VAF=sys.argv[3]	
    HF=sys.argv[4]
    MutAttMin=int(sys.argv[5])	
    MinVarError=5	
    ChaMat=Functions3.FillChaMat(ChaMat0,MutAttMin)
    print (ChaMat)
    Fas=ChaMat[:-4]+'_Position.fasta'	
    FinalRmLs=[Fas[:-6]+'.fa.json',Fas[:-6]+'.faOutSeq.fasta',Fas[:-6]+'.faOutSeqTopHap.fasta',Fas[:-6]+'_MutFre.txt',Fas[:-6]+'_Rm.txt',Fas[:-6]+'MOAIn.fasta',Fas[:-6]+'MOAlabelIn.txt']
    DataID=Fas[:-6]
    Functions3.ChaMat2Fas(ChaMat)
    if os.path.exists(Fas)!=True:	
        print ('making fasta file failed')
    else:
     CellLs,Cell2Seq=Functions3.ReadFasSeq(Fas)
     if float(VAF)>1:
            VAF=str(1.0*float(VAF)/len(CellLs))
     if float(HF)>1:
            HF=str(1.0*float(HF)/len(CellLs))
     if (float(HF)*len(CellLs))< MutAttMin:
         print ('use smaller number of mutation attach than ',(float(HF)*len(CellLs)))
     else:         
         print (VAF,HF,len(CellLs))
         Functions3.Fas2Json2(Fas)
         Json=Fas[:-3]+'.json'	
         if os.path.exists(Json)!=True: #count from 0
           print ('Json file was not created. Please check inout fasta file.')
         else:	
           OutSeq=Fas[:-3]+'OutSeq.fasta'	
           TopHapOut=Json[:-5]+'OutSeqTopHap.fasta'  
           Functions3.GetOut(TopHapOut,'>Outgroup\n'+open(OutSeq,'r').readlines()[0])   	   
           Hap=Json[:-5]+'Hap'	
           if os.path.exists(Hap)==True:
               print ('please delete ',Hap,' and try again')
           else:		   
               os.system(python+' vcf_json_parse.py '+Json+' --reference '+OutSeq+' --min_subgroup_size 3 --skip_mismatches --min_freq '+VAF+' -o '+Hap)		   
               Functions3.GetVAFpos(Hap+os.sep+'variant_stats.txt',float(VAF))	#Hap+os.sep+'variant_stats.txt'
               os.remove(Hap+os.sep+'Haplotypes.txt')
               os.remove(Hap+os.sep+'NA__2022_1.fasta')	
               os.remove(Hap+os.sep+'sequence_annotations.txt')
               os.remove(Hap+os.sep+'variant_stats.txt')		   
               os.system(python+' vcf_json_parse.py '+Json+' --reference '+OutSeq+' --min_subgroup_size 3 --skip_mismatches --min_freq 1.0 -o '+Hap+' --pos_file '+Hap+os.sep+'variant_stats_Force.txt') #position from 1			   
else:
    print ('input file is unclear. please select it by -Fas, -Json, or -Hap')
if os.path.exists(Hap)!=True:
    print ('haplotype alignment does not exists, so terminate the analysis')
else:
  VarSite=Hap+os.sep+'Haplotypes.txt'
  if len(open(VarSite,'r').readlines())>MaxSNV or len(open(VarSite,'r').readlines())<MinVarError:
    print ('Too many/little common sites were found', len(open(VarSite,'r').readlines()),'please change VAF and try again',VAF)
  else:	
    print ('do TopHap')
    os.system(python+' TopHap_main.py '+HF+' -Hap '+Hap+' -OutG '+TopHapOut)
    print (python+' TopHap_main.py '+HF+' -Hap '+Hap+' -OutG '+TopHapOut)
    OutFol=Hap+os.sep+'TopHapOut'+os.sep 
    if os.path.exists(OutFol)!=True:	
        os.mkdir(OutFol)	
    OutLs=glob.glob('TopHap*.csv')+glob.glob('TopHap*.nwk')+glob.glob('TopHap*.txt')+glob.glob('TopHap*.meg')+glob.glob('TopHap*.gv')
    OutLs+=['Time.txt','hf.txt','InputFasList.txt','TopHap.fasta']		
    for i in OutLs:
        shutil.copy(i, OutFol+i)
        os.remove(i)
    os.system('megacc -a analyze_user_tree_MP__nucleotide.mao -d '+OutFol+'TopHap.fasta -t '+OutFol+'TopHap_allMP_Anc.nwk -o '+OutFol+'TopHap_allMP_AncMP.nwk')		
    print ('map COI')
    print ('make input alignment')
    FileLs=open(OutFol+'InputFasList.txt','r').readlines()[1:]
    Out,OutSeq=Functions3.ReadFasSeq(Hap+os.sep+'OutG.fasta')	
    out=''#>Outgroup\n'+OutSeq[Out[0]]+'\n'
    Mout='ID\tsubregion\tcountry\tstate\tcollected\n'
    if sys.argv[1]!='-BEAM':	
      for HapIn in FileLs:
       HapLs,HapSeq=Functions3.ReadFasSeq(Hap+os.sep+HapIn.strip().replace('_FillAllCell',''))
       print ('File for COI computation',Hap+os.sep+HapIn.strip().replace('_FillAllCell',''))	   
       for H in HapLs:	   
         out+=H+'\n'+HapSeq[H]+'\n'
         Mout+=H[1:]+(('\t'+HapIn.strip())*3)+'\t20220412\n'
      MOAin=OutFol+'MOAfas.fasta'
      Functions3.GetOut(MOAin,out)	  
    if sys.argv[1]=='-BEAM':
      MOAin=BEAMinHap	
      OriFas=BEAMin
      CellLs0,Cell2Seq0=Functions3.ReadMegSeq(OriFas)
      RefFas=Meg
    else: 
      OriFas=Fas
      CellLs0,Cell2Seq0=Functions3.ReadFasSeq(OriFas)
      RefFas=Fas      
    SNVc=len(Cell2Seq0[CellLs0[0]])  
    COImatrix=OutFol+'MOAout'+os.sep+'COI_matrix.txt'   
    GVin=COImatrix.replace('MOAout'+os.sep+'COI_matrix.txt','TopHap_MutTree.gv')   
    Functions3.MakeMOAin1(Hap+os.sep+'Haplotypes.txt',GVin,MOAin)
    MOAin1=MOAin[:-6]+'MOAIn.fasta'	
    print ('compute COI')
    os.system(python+' FastMOA_update.py '+MOAin1+' '+Hap+os.sep+os.sep+'HaplotypesMOAlabelIn.txt  -o '+OutFol+'MOAout --threshold 0.1 --disable_graph_flipping --flip_pass_thresh 1.1')
    print ('map COI')
    print (MOAin1)
    os.system(python+ ' MapCOITest.py '+COImatrix+' '+MOAin1+' '+Hap+os.sep+'Haplotypes.txt '+OriFas+ ' '+str(COIcut))
    print (python+ ' MapCOITest.py '+COImatrix+' '+MOAin1+' '+Hap+os.sep+'Haplotypes.txt '+OriFas+ ' '+str(COIcut))# '+' '+python+' '+MLTED)	
    print (MOAin1)
    FPFN=open(OutFol+'FPFN.txt','r').readlines()
    FP=float(FPFN[0].split(':')[1])
    FN=float(FPFN[1].split(':')[1])
    print ('FP FN',FP,FN)
    if sys.argv.count('-MutID')!=0: 
       MutIDFile=sys.argv[sys.argv.index('-MutID')+1]
    else: MutIDFile='NA'
    print (MutIDFile,sys.argv.count('-MutID'),sys.argv)   
    MutIDFile='NA'    
    Cont='y'
    while Cont=='y':
        print ('final COI test for tip recurrent and backward')
        if os.path.exists(OutFol+'MutTree_Final.gv')!=True: IniGV=OutFol+'TopHap_MutTree_prune_COI3_ave3_rec_back.gv'#OutFol+'TopHapCOI_MutTree_COI3_ave3.gv'
        else: IniGV=OutFol+'MutTree_Final.gv'
        Functions3.PruneTip(IniGV,OutFol+'TopHap.fasta',OutFol+'TopHap_allMP_rooted.nwk',OutFol+'TopHapNA__2022_1_FillAllCell_anno1.txt')
        Functions3.CleanMutTree(IniGV[:-3]+'1.gv','NA',OutFol+'TopHap_allMP_rooted.nwk') #output 'MutTree.gv'      
        Cont=Functions3.AnaTipRec(OutFol+'MutTree.gv',COImatrix,OriFas,FP,FN) #output 'MutTree_Final.gv'
        print (Cont)
    print  ('Final')# CleanMutTree.py '+OutFol+'MutTree_Final.gv '+MutIDFile)   
    Functions3.CleanMutTree(OutFol+'MutTree_Final.gv',MutIDFile,OutFol+'TopHap_allMP_rooted.nwk')    #'MutTree.gv'  
    print ('MOA')
    os.system(python+' MapMissMut.py '+OutFol+' '+OutFol+'MutTree.gv '+OriFas+' '+RefFas+' '+str(MOAcut)+' 3 3 python') #3 2
    print (python+' MapMissMut.py '+OutFol+' '+OutFol+'MutTree.gv '+OriFas+' '+RefFas+' '+str(MOAcut)+' 3 3 python')
    print ('map back recurrent and back mutations')
    TopHapMOAgv= OutFol+'MOAprune_TopHap.gv'# glob.glob(OutFol+'MOAprune_*TopHap0_TopHap.gv')[0] 
    TopHapMOACellAnno= OutFol+'MOAprune_TopHap_FilInter_CellAnnoAll_final_CellAnnoAll_final.txt'#glob.glob(OutFol+'MOAprune_*TopHap0_TopHap_FilInter_CellAnnoAll_final_CellAnnoAll_final.txt')[0]
    TopHapMOAfas=OutFol+'MOAprune_TopHap_FilInter_CellAnnoAll_final_CellAnnoAll_final.fasta'#glob.glob(OutFol+'MOAprune_TopHap0*_TopHap_FilInter_CellAnnoAll_final_CellAnnoAll_final.fasta')[0]       
    Functions3.TestBackRec2(OutFol+'MutTree.gv',TopHapMOAgv,TopHapMOACellAnno,OriFas,FP,FN)  
    print ('Functions in',OutFol+'MutTree.gv',TopHapMOAgv,TopHapMOACellAnno,OriFas,FP,FN)    
    TopHapMOARecBackgv=TopHapMOAgv[:-3]+'_backrec.gv'
    if sys.argv.count('-MutID')!=0: 
       MutIDFile=sys.argv[sys.argv.index('-MutID')+1]
    else: MutIDFile='NA'
    print ('finalize mutation tree and map mutation ID')    
    Functions3.CleanMutTree1(TopHapMOARecBackgv,MutIDFile,TopHapMOACellAnno,SNVc) #output '_final.gv'     
    if os.path.exists(OutFol+'TopHapCOIMOA0.nwk')!=True: os.system('megacc -a infer_MP_nucleotide.mao -d ' + TopHapMOAfas +' -o ' + OutFol+'TopHapCOIMOA0.nwk')   
    if os.path.exists(OutFol+'TopHapCOIMOA.nwk')!=True: os.system('megacc -a analyze_user_tree_MP__nucleotide.mao -d '+TopHapMOAfas+' -t '+OutFol+'TopHapCOIMOA0.nwk -o '+OutFol+'TopHapCOIMOA.nwk')	
    Keep=[OutFol+'FPFN.txt',OutFol+'TopHap_allMP_rooted_1_prune.nwk',OutFol+'MutTree.png',OutFol+'MutTree.gv']
    Keep+=[TopHapMOAfas,TopHapMOACellAnno,TopHapMOARecBackgv[:-3]+'_final.gv',OutFol+'TopHapCOIMOA.nwk']
    Keep+=[OutFol+'TopHapNA__2022_1_FillAllCell_anno1.txt',OutFol+'TopHap.fasta',OutFol+'TopHap_allMP_AncMP.nwk',OutFol+'TopHap_MutTree.gv']#,OutFol+'MutTree_COI.png'] #,OutFol+'TopHap_MutTree_COI3_ave3.png',OutFol+'TopHap_MutTree_COI3_ave3.gv'
    Keep+=[OutFol+'COIsummary1.txt',OutFol+'FPFN.txt',OutFol+'RecBackSupp.txt']	
    Keep+=[OutFol+'MOAout'+os.sep+'COI_matrix.txt',OutFol+'MOAout'+os.sep+'COI.txt',OutFol+'MOAout'+os.sep+'co_counts_matrix.txt']	#OutFol+'MOAout'+os.sep+'COI_matrix.txt',
    AllFile=glob.glob(OutFol+'*.txt')+glob.glob(OutFol+'*.gv')+glob.glob(OutFol+'*.fasta')+glob.glob(OutFol+'*.nwk')+glob.glob(OutFol+'*.meg')+glob.glob(OutFol+'*.csv')
    MOAOutFol=glob.glob(OutFol+'MOAmap*')
    for mFol in MOAOutFol:
        Clean(mFol+os.sep+'*.fas')
        Clean(mFol+os.sep+'*.txt')
        Clean(mFol+os.sep+'*.png')
        os.rmdir(mFol)    
    AllFile=list(set(AllFile))
    for CFile in AllFile:
       if Keep.count(CFile)==0 : os.remove(CFile)
    Clean(DataID+'_BEAM_HapMOAIn.fasta')  
    Clean(DataID+'_BEAM.fa.json')   
    Clean(DataID+'_BEAM.faOutSeq.fasta') 
    Clean(DataID+'_BEAM.faOutSeqTopHap.fasta')
    Clean(Hap+os.sep+'HaplotypesMOAlabelIn.txt')
    Clean(Hap+os.sep+'*.fasta')
    Clean(Hap+os.sep+'variant_*.txt')  
    CopyAndDel(OutFol+'TopHap_allMP_rooted_1_prune.nwk',OutFol+'TopHapCOI.nwk')
    CopyAndDel(OutFol+'MutTree.png',OutFol+'TopHapCOI.png')
    CopyAndDel(OutFol+'MutTree.gv',OutFol+'TopHapCOI.gv')
    THfas=OutFol+'TopHapCOIMOA.fasta'
    THclo=OutFol+'TopHapCOIMOA.txt'
    CopyAndDel(TopHapMOAfas,OutFol+'TopHapCOIMOA.fasta')
    CopyAndDel(TopHapMOACellAnno,OutFol+'TopHapCOIMOA.txt')
    CopyAndDel(TopHapMOARecBackgv[:-3]+'_final.gv',OutFol+'TopHapCOIMOA.gv')
    CopyAndDel(OutFol+'TopHapNA__2022_1_FillAllCell_anno1.txt',OutFol+'TopHap.txt')
    if os.path.exists(OutFol+'TopHap_allMP_AncMP.nwk')==True: CopyAndDel(OutFol+'TopHap_allMP_AncMP.nwk',OutFol+'TopHap.nwk')
    elif os.path.exists(OutFol+'TopHap_allMP_Anc.nwk')==True: CopyAndDel(OutFol+'TopHap_allMP_Anc.nwk',OutFol+'TopHap.nwk')
    CopyAndDel(OutFol+'TopHap_MutTree.gv',OutFol+'TopHap.gv')
  #  os.system('dot -Tpng '+OutFol+'TopHap.gv'+' -o '+OutFol+'TopHap.png') 
    os.system('dot -Tpng '+OutFol+'TopHapCOIMOA.gv'+' -o '+OutFol+'TopHapCOIMOA.png')     
    Functions3.GetOut(OutFol+'Param.txt','\n'.join(['COIcut: '+str(COIcut),'MOAcut: '+str(MOAcut),'VAF: '+str(VAF),'HF: '+str(HF)]))
    CopyAndDel(OutFol+'TopHapCOIMOA.gv',Hap+os.sep+'TopHapPlus.gv')   
    CopyAndDel(OutFol+'TopHapCOIMOA.png',Hap+os.sep+'TopHapPlus.png')   
    CopyAndDel(OutFol+'TopHapCOIMOA.txt',Hap+os.sep+'TopHapPlus.txt')    
   # CopyAndDel(OutFol+'TopHapCOIMOA.fasta',Hap+os.sep+'TopHapPlus.fasta')
    CopyAndDel(OutFol+'TopHapCOIMOA.nwk',Hap+os.sep+'TopHapPlus.nwk')  
    Clean(OutFol+'MOAout'+os.sep+'*') 
    os.rmdir(OutFol+'MOAout')     
    Clean(OutFol+'*')  
    os.rmdir(OutFol)     
    AllFileLs=glob.glob(Hap+os.sep+'*') #+glob.glob(OutFol+'*.gv')+glob.glob(OutFol+'*.png')+glob.glob(OutFol+'*.fasta')+glob.glob(OutFol+'*.nwk')
    for i in AllFileLs:
         if i.find(Hap+os.sep+'TopHapPlus.')==-1: os.remove(i)  
    for i in FinalRmLs:
         if os.path.exists(i)==True: os.remove(i)    
    if sys.argv.count('-ST')!=0:
        print ('all clone for ST')
        MutPls=BEAMin[:-4]+'MOAlabelIn.txt'
        Functions3.SeqAnnoST(THfas,THclo,AllCellFas,MutPls)     