 !*** INPUT generated by FMOPhore, Copyright "©" reserved to Peter Ezzat Girgis Fahmy Ibrahim. ***
 $CONTRL RUNTYP=ENERGY NPRINT=-5 ISPHER=-1 MAXIT=50 $END
 $SYSTEM MWORDS=29 $END
 $GDDI NGROUP=1 $END
 $INTGRL NINTIC=-9000000 $END
 $SCF DIRSCF=.TRUE. NPUNCH=0 $END
 $BASIS GBASIS=DFTB $END
 $PCM SOLVNT=WATER ICOMP=2 ICAV=1 IDISP=1 IFMO=-1 $END
 $PCMCAV RADII=SUAHF $END
 $TESCAV NTSALL=60 $END
 $DFTB
   SCC=.TRUE. DFTB3=.TRUE. DAMPXH=.TRUE. DAMPEX=4.00 DISP=UFF
 $END
 $FMOPRP
    NAODIR=200
    NGRFMO(1)=1, 1, 0, 0, 0,   0, 0, 0, 0, 0
    IPIEDA=1
    NPRINT=9
    modorb=3
 $END
 $DFTBSK
   Br Br /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Br-Br.skf
   Br C /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Br-C.skf
   Br Ca /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Br-Ca.skf
   Br Cl /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Br-Cl.skf
   Br F /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Br-F.skf
   Br H /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Br-H.skf
   Br I /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Br-I.skf
   Br K /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Br-K.skf
   Br Mg /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Br-Mg.skf
   Br N /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Br-N.skf
   Br Na /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Br-Na.skf
   Br O /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Br-O.skf
   Br P /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Br-P.skf
   Br S /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Br-S.skf
   Br Zn /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Br-Zn.skf
   C Br /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/C-Br.skf
   C C /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/C-C.skf
   C Ca /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/C-Ca.skf
   C Cl /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/C-Cl.skf
   C F /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/C-F.skf
   C H /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/C-H.skf
   C I /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/C-I.skf
   C K /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/C-K.skf
   C Mg /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/C-Mg.skf
   C N /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/C-N.skf
   C Na /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/C-Na.skf
   C O /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/C-O.skf
   C P /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/C-P.skf
   C S /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/C-S.skf
   C Zn /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/C-Zn.skf
   Ca Br /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Ca-Br.skf
   Ca C /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Ca-C.skf
   Ca Ca /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Ca-Ca.skf
   Ca Cl /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Ca-Cl.skf
   Ca F /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Ca-F.skf
   Ca H /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Ca-H.skf
   Ca I /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Ca-I.skf
   Ca K /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Ca-K.skf
   Ca Mg /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Ca-Mg.skf
   Ca N /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Ca-N.skf
   Ca Na /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Ca-Na.skf
   Ca O /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Ca-O.skf
   Ca P /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Ca-P.skf
   Ca S /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Ca-S.skf
   Ca Zn /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Ca-Zn.skf
   Cl Br /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Cl-Br.skf
   Cl C /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Cl-C.skf
   Cl Ca /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Cl-Ca.skf
   Cl Cl /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Cl-Cl.skf
   Cl F /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Cl-F.skf
   Cl H /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Cl-H.skf
   Cl I /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Cl-I.skf
   Cl K /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Cl-K.skf
   Cl Mg /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Cl-Mg.skf
   Cl N /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Cl-N.skf
   Cl Na /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Cl-Na.skf
   Cl O /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Cl-O.skf
   Cl P /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Cl-P.skf
   Cl S /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Cl-S.skf
   Cl Zn /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Cl-Zn.skf
   F Br /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/F-Br.skf
   F C /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/F-C.skf
   F Ca /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/F-Ca.skf
   F Cl /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/F-Cl.skf
   F F /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/F-F.skf
   F H /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/F-H.skf
   F I /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/F-I.skf
   F K /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/F-K.skf
   F Mg /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/F-Mg.skf
   F N /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/F-N.skf
   F Na /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/F-Na.skf
   F O /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/F-O.skf
   F P /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/F-P.skf
   F S /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/F-S.skf
   F Zn /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/F-Zn.skf
   H Br /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/H-Br.skf
   H C /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/H-C.skf
   H Ca /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/H-Ca.skf
   H Cl /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/H-Cl.skf
   H F /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/H-F.skf
   H H /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/H-H.skf
   H I /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/H-I.skf
   H K /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/H-K.skf
   H Mg /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/H-Mg.skf
   H N /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/H-N.skf
   H Na /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/H-Na.skf
   H O /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/H-O.skf
   H P /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/H-P.skf
   H S /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/H-S.skf
   H Zn /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/H-Zn.skf
   I Br /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/I-Br.skf
   I C /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/I-C.skf
   I Ca /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/I-Ca.skf
   I Cl /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/I-Cl.skf
   I F /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/I-F.skf
   I H /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/I-H.skf
   I I /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/I-I.skf
   I K /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/I-K.skf
   I Mg /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/I-Mg.skf
   I N /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/I-N.skf
   I Na /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/I-Na.skf
   I O /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/I-O.skf
   I P /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/I-P.skf
   I S /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/I-S.skf
   I Zn /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/I-Zn.skf
   K Br /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/K-Br.skf
   K C /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/K-C.skf
   K Ca /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/K-Ca.skf
   K Cl /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/K-Cl.skf
   K F /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/K-F.skf
   K H /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/K-H.skf
   K I /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/K-I.skf
   K K /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/K-K.skf
   K Mg /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/K-Mg.skf
   K N /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/K-N.skf
   K Na /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/K-Na.skf
   K O /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/K-O.skf
   K P /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/K-P.skf
   K S /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/K-S.skf
   K Zn /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/K-Zn.skf
   Mg Br /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Mg-Br.skf
   Mg C /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Mg-C.skf
   Mg Ca /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Mg-Ca.skf
   Mg Cl /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Mg-Cl.skf
   Mg F /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Mg-F.skf
   Mg H /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Mg-H.skf
   Mg I /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Mg-I.skf
   Mg K /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Mg-K.skf
   Mg Mg /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Mg-Mg.skf
   Mg N /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Mg-N.skf
   Mg Na /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Mg-Na.skf
   Mg O /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Mg-O.skf
   Mg P /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Mg-P.skf
   Mg S /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Mg-S.skf
   Mg Zn /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Mg-Zn.skf
   N Br /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/N-Br.skf
   N C /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/N-C.skf
   N Ca /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/N-Ca.skf
   N Cl /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/N-Cl.skf
   N F /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/N-F.skf
   N H /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/N-H.skf
   N I /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/N-I.skf
   N K /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/N-K.skf
   N Mg /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/N-Mg.skf
   N N /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/N-N.skf
   N Na /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/N-Na.skf
   N O /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/N-O.skf
   N P /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/N-P.skf
   N S /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/N-S.skf
   N Zn /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/N-Zn.skf
   Na Br /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Na-Br.skf
   Na C /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Na-C.skf
   Na Ca /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Na-Ca.skf
   Na Cl /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Na-Cl.skf
   Na F /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Na-F.skf
   Na H /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Na-H.skf
   Na I /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Na-I.skf
   Na K /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Na-K.skf
   Na Mg /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Na-Mg.skf
   Na N /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Na-N.skf
   Na Na /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Na-Na.skf
   Na O /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Na-O.skf
   Na P /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Na-P.skf
   Na S /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Na-S.skf
   Na Zn /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Na-Zn.skf
   O Br /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/O-Br.skf
   O C /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/O-C.skf
   O Ca /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/O-Ca.skf
   O Cl /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/O-Cl.skf
   O F /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/O-F.skf
   O H /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/O-H.skf
   O I /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/O-I.skf
   O K /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/O-K.skf
   O Mg /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/O-Mg.skf
   O N /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/O-N.skf
   O Na /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/O-Na.skf
   O O /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/O-O.skf
   O P /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/O-P.skf
   O S /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/O-S.skf
   O Zn /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/O-Zn.skf
   P Br /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/P-Br.skf
   P C /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/P-C.skf
   P Ca /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/P-Ca.skf
   P Cl /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/P-Cl.skf
   P F /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/P-F.skf
   P H /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/P-H.skf
   P I /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/P-I.skf
   P K /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/P-K.skf
   P Mg /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/P-Mg.skf
   P N /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/P-N.skf
   P Na /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/P-Na.skf
   P O /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/P-O.skf
   P P /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/P-P.skf
   P S /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/P-S.skf
   P Zn /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/P-Zn.skf
   S Br /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/S-Br.skf
   S C /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/S-C.skf
   S Ca /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/S-Ca.skf
   S Cl /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/S-Cl.skf
   S F /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/S-F.skf
   S H /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/S-H.skf
   S I /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/S-I.skf
   S K /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/S-K.skf
   S Mg /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/S-Mg.skf
   S N /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/S-N.skf
   S Na /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/S-Na.skf
   S O /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/S-O.skf
   S P /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/S-P.skf
   S S /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/S-S.skf
   S Zn /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/S-Zn.skf
   Zn Br /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Zn-Br.skf
   Zn C /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Zn-C.skf
   Zn Ca /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Zn-Ca.skf
   Zn Cl /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Zn-Cl.skf
   Zn F /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Zn-F.skf
   Zn H /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Zn-H.skf
   Zn I /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Zn-I.skf
   Zn K /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Zn-K.skf
   Zn Mg /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Zn-Mg.skf
   Zn N /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Zn-N.skf
   Zn Na /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Zn-Na.skf
   Zn O /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Zn-O.skf
   Zn P /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Zn-P.skf
   Zn S /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Zn-S.skf
   Zn Zn /PATH-TO-GAMES-SOFTWARE/software/gamess/auxdata/DFTB/3ob-3-1/Zn-Zn.skf
 $END
 $FMO
      SCFTYP(1)=RHF
      MODGRD=10
      MODMUL=0
      respap=0
      resppc=4
      NBODY=2
      MAXCAO=5

 $END
 $FMOLMO
 DFTB-C 4 4
   1 0  0.558756  0.000000  0.000000  0.829332
   0 1  0.558757  0.781901  0.000000 -0.276445
   0 1  0.558756 -0.390951  0.677146 -0.276445
   0 1  0.558756 -0.390951 -0.677146 -0.276445
 $END
 $FMOBND

 C1 
 C    6
 O    8
 H    1
 N    7
 S   16
 F    9
 P   15
 Cl  17
 Br  35
 I   53
 $END
 $FMOXYZ