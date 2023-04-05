import Functions3
import os
import glob
import sys
import time

MinShare=10 #obs>= :include to TopHap if hf is not useful
StartTime='Start time: '+time.strftime("%Y-%m-%d %H:%M")+'\n'
THCutFre=float(sys.argv[1])#0.05 #hf
HFabsMinMax=[1,1000]

MAO='infer_ML_nucleotide.mao'#'infer_NJ_nucleotide.mao'#'infer_ML_nucleotide.mao'#'infer_MP_nucleotide.mao'

Python='python'
OutGfas='Outgroup.fasta'
AnnoFile=''#'C:\\Users\\tuf78332\\Desktop\\TopHap\\FRACTAL\\Hap0.25_Top50\\Hap0.25_Top50\\TopHapNA_2022_1_anno.txt'
PosStart=0
path = os.path.join(os.path.expanduser('~'), os.getcwd())
Fill='y'

if os.path.exists('Bootstrap')==True:
     print ('Please delete Bootstrap directory and try again')
      

else:

 if sys.argv.count('-Hap')==1:
    Fol=sys.argv[sys.argv.index('-Hap')+1]+os.sep
    FasLs00=glob.glob(Fol+'*.fasta')
    out='Input\n'	
    for Fas in FasLs00:
       print ('input: ',Fas)	
       InLs,InSeq=Functions3.ReadFasSeq(Fas)	
       if Fas!=Fol+'OutG.fasta' and len(InLs)>=1:	
          if Fill=='y': out+=Fas.split(os.sep)[-1][:-6]+'_FillAllCell.fasta'+'\n'
          else: out+=Fas.split(os.sep)[-1]+'\n'		  
       if Fill=='y':
          if os.path.exists(Fas.split(os.sep)[-1][:-6]+'_FillAllCell.fasta')!=True: 
              Functions3.FillMissBase_consen(Fas, 'A,T,C,G',0.9)#os.system(Python+' FillMissBase_consen.py '+Fas+ ' A,T,C,G 0.9')
              os.remove('UniqueSeq.fasta')           			   

    if Fill=='y': FasLs0=glob.glob(Fol+'*_FillAllCell.fasta')	
    else: FasLs0=glob.glob(Fol+'*.fasta')		
    Functions3.GetOut('InputFasList.txt',out)
    FasTaLs=['InputFasList.txt']

 if sys.argv.count('-OutG')==1: 
     OutGfas=sys.argv[sys.argv.index('-OutG')+1]
 else: print ('please provide the directory of haplotype alignments')	

 PosLs0=open(Fol+'Haplotypes.txt','r').readlines()
 PosLs=[]
 for i in PosLs0: 
       if PosStart==0: PosLs.append(int(i)+1)
       else: PosLs.append(int(i))	
 Ols,O2Seq=Functions3.ReadFasSeq(OutGfas)
 outG=''
 OseqLs=[]
 for O in Ols:
    Seq=O2Seq[O]
    Hap=''	
    for Pos in PosLs:
      	
       Hap+=Seq[Pos-1]
    outG +=O+'\n'+Hap+'\n'
    OseqLs.append(Hap)	


 TopHap2Data={}
 Sharedcount=[]
 for FasLs in FasTaLs:
    print ('input: ',FasLs)
    Name=''	
    FasLs=open(FasLs,'r').readlines()[1:]
    for i in FasLs:
            print ('making TopHap sequences for the slice',i)	
            Fas0=i.strip().split('\t')[0] ####Fix
            Fas=Fol+Fas0	 ####Fix		
            Seq2Cou={}
            StLs,St2Seq=Functions3.ReadFasSeq(Fas)
            Fill=1.0*len(StLs)*THCutFre			
            for i in StLs:
                       Seq2Cou[St2Seq[i]]=Seq2Cou.get(St2Seq[i],0)+1	
            CouLs=Functions3.InvertDic2Ls(Seq2Cou)
           
            CouLs.sort(reverse=True)
            print (CouLs,Seq2Cou)
          
            for Seq in Seq2Cou:
					 
                  if Seq2Cou[Seq]>=Fill :
                    if OseqLs.count(Seq)==0:
                       				
                       TopHap2Data[Seq]=TopHap2Data.get(Seq,[])+[Fas0]	####Fix	
                   					   
            THc=len(TopHap2Data)
            c=0	
            Cou=CouLs[c]	
            			
            while THc<HFabsMinMax[0] and Cou>=MinShare and c<len(CouLs) :
               Cou=CouLs[c]
             
               if Cou>=MinShare:
                 Sharedcount.append(Fas0+' hf reduced to abs_hf '+str(Cou))				   
                 for Seq in Seq2Cou:
                  if Seq2Cou[Seq]==Cou:				 

                    if Seq not in TopHap2Data:
                      if OseqLs.count(Seq)==0:
				  
                          TopHap2Data[Seq]=TopHap2Data.get(Seq,[])+[Fas0]	####Fix	
                    						  
                 THc=len(TopHap2Data)
                				 
               c+=1				 
 print ('tophap c',len(TopHap2Data))
 Functions3.GetOut('hf.txt','\n'.join(Sharedcount)+'\nMinimum number of TopHap haplotypes to extract: '+str(HFabsMinMax[0]))
 

 if len(TopHap2Data)>HFabsMinMax[1]:
    print ('>',HFabsMinMax[1],' sequences are found, please use smaller value for HFabsMinMax')
    print (CouLs)	

 elif len(TopHap2Data)<3:
    print (len(TopHap2Data),' no TopHap sequences produced')
    print (CouLs)	
 else:
  TopID=1
  outF=''
  ID2fasID={}
  TopHapLs=list(TopHap2Data.keys())
  FasTopHap0='TopHap_withMiss.fasta'
  out=''
  c=1
  for i in TopHapLs:
     out+='>'+str(c)+'\n'+i+'\n'
     c+=1	
  Functions3.GetOut(FasTopHap0,out)	 
  Functions3.FillMissBase_consen(FasTopHap0, 'A,T,C,G',1.0)#os.system(Python+' FillMissBase_consen.py '+FasTopHap0+ ' A,T,C,G 1.0')
  AA,TopHap02Seq=Functions3.ReadFasSeq(FasTopHap0[:-6]+'_FillAllCell.fasta')
  UniqueMissCellLs,MisscellSeq=Functions3.ReadFasSeq('UniqueSeq.fasta')
  TopHapSeqLs=[]
  for i in TopHap02Seq:
   if TopHap02Seq[i].find('?')==-1: 
    TopHapSeqLs.append(TopHap02Seq[i])
  for i in MisscellSeq:
    print ('seq with miss',i)
    TopHapSeqLs.append(MisscellSeq[i])	
  TopHapSeqLs=list(set(TopHapSeqLs))
  os.remove('UniqueSeq.fasta')
  os.remove(FasTopHap0)
  os.remove(FasTopHap0[:-6]+'_FillAllCell.fasta')

  for Hap in TopHapSeqLs:
    ID='All'+str(TopID)	
    if OseqLs.count(Hap)==0:	
      outF+='>'+ID+'\n'+Hap+'\n'
      TopID+=1
  outF=outG+outF	
  Functions3.GetOut(Fol+'OutG.fasta',outG)	
  Functions3.GetOut('TopHap.fasta',outF)

  if os.path.exists(AnnoFile)!=True:
   print ('annotating haplotypes...')
   for Fas in FasLs0:
     Functions3.TopHapAnnotation('TopHap.fasta',Fas)
     Functions3.AnnotateOriID('TopHap'+Fas[:-6].split(os.sep)[-1]+'_anno.txt',Fol+'sequence_annotations.txt')#os.system(Python+' AnnotateOriID.py TopHap'+Fas[:-6].split(os.sep)[-1]+'_anno.txt '+Fol+'sequence_annotations.txt')
     if os.path.exists('TopHap'+Fas[:-6].split(os.sep)[-1]+'_anno1.txt')!=True:
         print ('AnnotateOriID.py failed')
	 
      	
  print ('infer tree')
  os.system('megacc -a '+MAO+' -d ' + 'TopHap.fasta' +' -o ' + 'TopHap_allMP.nwk')
  AncFileLs=glob.glob('TopHap_allMP_ancestral_states_*.txt')
  for AncF in AncFileLs:
               os.remove(AncF)	  
  if os.path.exists('TopHap_allMP_summary.txt')==True: os.remove('TopHap_allMP_summary.txt')	
  print ('infer mutation tree')
  Functions3.TopHap_DrawMutTree2('TopHap.fasta',Fol+'Haplotypes.txt','TopHap_allMP.nwk','Outgroup')
 # os.system(Python+' TopHap_DrawMutTree2.py TopHap.fasta '+Fol+'Haplotypes.txt TopHap_allMP.nwk Outgroup '+Python)
  if os.path.exists('TopHap_allMP.nwk')!=True:
    print ('tree was not produced. See MEGA error message')



  os.remove('TopHap_allMP.nwk') 
EndTime='End time: '+time.strftime("%Y-%m-%d %H:%M")+'\n'
Functions3.GetOut('Time.txt',StartTime+EndTime)	
