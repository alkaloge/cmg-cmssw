#!/bin/bash

if [[ "$1" == "afs" ]]; then
    T="/afs/cern.ch/work/g/gpetrucc/TREES_72X_040115";
    J=4;
elif [[ "$1" == "SingleLepAFS" ]]; then
#    T="/afs/cern.ch/work/k/kirschen/public/PlotExampleSamples/V2";
    T="/afs/cern.ch/work/k/kirschen/public/PlotExampleSamples/V3";
    J=4;
elif [[ "$HOSTNAME" == "cmsphys10" ]]; then
    T="/data/g/gpetrucc/TREES_72X_040115";
    J=8;
else
    T="/afs/cern.ch/work/g/gpetrucc/TREES_72X_040115";
    J=4;
fi

LUMI=4.0
OUTDIR="susy_cards_1l_4fb"
#OPTIONS=" -P $T -j $J -l $LUMI -f --s2v --tree treeProducerSusyMultilepton --od $OUTDIR --asimov "
#OPTIONS=" $OPTIONS -F sf/t $T/0_lepMVA_v1/evVarFriend_{cname}.root "
OPTIONS=" -P $T -j $J -l $LUMI -f --s2v --tree treeProducerSusySingleLepton --od $OUTDIR --asimov "


function makeCard_2lss {
    local EXPR=$1; local BINS=$2; local SYSTS=$3; local OUT=$4; local GO=$5

    # b-jet cuts
    case $SR in
    0[0-9X])  GO="${GO} -R nBjet nBjet0 nBJetMedium40==0 " ;;
    1[0-9X])  GO="${GO} -R nBjet nBjet1 nBJetMedium40==1 " ;;
    2[0-9X])  GO="${GO} -R nBjet nBjet2 nBJetMedium40==2 " ;;
    3[0-9X])  GO="${GO} -R nBjet nBjet3 nBJetMedium40>=3 " ;;
    2[0-9X]+)  GO="${GO} -R nBjet nBjet2 nBJetMedium40>=2 " ;;
    0[0-9X]s)  GO="${GO} -R nBjet nBjet3 nBJetMedium40+min(nBJetMedium25+nSoftBTight25-nBJetMedium40,1)==0 " ;;
    1[0-9X]s)  GO="${GO} -R nBjet nBjet3 nBJetMedium40+min(nBJetMedium25+nSoftBTight25-nBJetMedium40,1)==1 " ;;
    2[0-9X]s)  GO="${GO} -R nBjet nBjet3 nBJetMedium40+min(nBJetMedium25+nSoftBTight25-nBJetMedium40,1)==2 " ;;
    3[0-9X]s)  GO="${GO} -R nBjet nBjet3 nBJetMedium40+min(nBJetMedium25+nSoftBTight25-nBJetMedium40,1)>=3 " ;;
    esac;

    # kinematics
    case $SR in
    [0-3]X|[0-3]X+)  GO="${GO} -R met met met_pt>50 -R ht ht htJet40j>200 " ;;
    esac;

    # lepton final state
    case $LL in
    ee)  GO="${GO} -R anyll ee abs(LepGood1_pdgId)==11&&abs(LepGood2_pdgId)==11 " ;;
    em)  GO="${GO} -R anyll em abs(LepGood1_pdgId)!=abs(LepGood2_pdgId) " ;;
    mm)  GO="${GO} -R anyll mm abs(LepGood1_pdgId)==13&&abs(LepGood2_pdgId)==13 " ;;
    3l)  GO="${GO} -I exclusive -X same-sign -R anyll lep3-cuts LepGood3_relIso03<0.1&&LepGood3_tightId>(abs(LepGood3_pdgId)==11)&&LepGood3_sip3d<4&&(abs(LepGood3_pdgId)==13||(LepGood3_convVeto&&LepGood3_lostHits==0&&LepGood3_tightCharge>1))"
    esac;

    # lepton pt categories
    case $LPt in
    hl)  GO="${GO} -I lep2_pt25" ;;
    ll)  GO="${GO} -I lep1_pt25 -X lep2_pt25" ;;
    ii)  GO="${GO} -X lep1_pt25 -X lep2_pt25" ;;
    2020)  GO="${GO} -R lep1_pt25 lep2020 LepGood2_pt>20 -X lep2_pt25" ;;
    esac;

    # inclusive vs exclusive
    case $MOD in
    inc) GO="${GO} -X exclusive --mcc bins/susymultilepton/susy_2lssinc_lepchoice.txt" ;;
    esac;

    if [[ "$PRETEND" == "1" ]]; then
        echo "making datacard $OUT from makeShapeCardsSusy.py mca-Phys14.txt bins/susymultilepton/susy_2lss_sync.txt \"$EXPR\" \"$BINS\" $SYSTS $GO;"
    else
        echo "making datacard $OUT from makeShapeCardsSusy.py mca-Phys14.txt bins/susymultilepton/susy_2lss_sync.txt \"$EXPR\" \"$BINS\" $SYSTS $GO;"
        python makeShapeCardsSusy.py mca-Phys14.txt bins/susymultilepton/susy_2lss_sync.txt "$EXPR" "$BINS" $SYSTS -o $OUT $GO;
        echo "  -- done at $(date)";
    fi;
}


function makeCard_1l {
    local EXPR=$1; local BINS=$2; local SYSTS=$3; local OUT=$4; local GO=$5

    # b-jet cuts
    case $nB in
    0B)  GO="${GO} -R 1nB 0nB nBJetCMVAMedium30==0 " ;;
    1B)  GO="${GO} -R 1nB 1nB nBJetCMVAMedium30==1 " ;;
    2B)  GO="${GO} -R 1nB 2nB nBJetCMVAMedium30==2 " ;;
    2Btop)  GO="${GO} -R 1nB 2nB nBJetCMVAMedium30==2&&Topness>5 " ;;
    3p)  GO="${GO} -R 1nB 3nBp nBJetCMVAMedium30>=3 " ;;
    esac;

    # lepton final state
    case $ST in
    ST0)  GO="${GO} -R st200 st200250 LepGood_pt[tightLeps25idx[0]]+met_pt>200&&LepGood_pt[tightLeps25idx[0]]+met_pt<250 " ;;
    ST1)  GO="${GO} -R st200 st250350 LepGood_pt[tightLeps25idx[0]]+met_pt>250&&LepGood_pt[tightLeps25idx[0]]+met_pt<350 " ;;
    ST2)  GO="${GO} -R st200 st350450 LepGood_pt[tightLeps25idx[0]]+met_pt>350&&LepGood_pt[tightLeps25idx[0]]+met_pt<450 " ;;
    ST3)  GO="${GO} -R st200 st450600 LepGood_pt[tightLeps25idx[0]]+met_pt>450&&LepGood_pt[tightLeps25idx[0]]+met_pt<600 " ;;
    ST4)  GO="${GO} -R st200 st600Inf LepGood_pt[tightLeps25idx[0]]+met_pt>600 " ;;
    esac;

    # lepton pt categories
    case $nJ in
    45j)  GO="${GO} -R geq6j 45j nCentralJet30>=4&&nCentralJet30<=5"  ;;
    68j)  GO="${GO} -R geq6j 67j nCentralJet30>=6&&nCentralJet30<=8"  ;;
    6Infj)  GO="${GO} -R geq6j geq6j nCentralJet30>=6"  ;;
    9Infj)  GO="${GO} -R geq6j geq8j nCentralJet30>=9"  ;;
    68TTj)  GO="${GO} -R geq6j 79TTj nCentralJet30+2*nHighPtTopTagPlusTau23>=6&&nCentralJet30+2*nHighPtTopTagPlusTau23<9"  ;; ##-R dphi --> 0.5 und SingleTopness
    9InfTTj)  GO="${GO} -R geq6j 9InfTTj nCentralJet30+2*nHighPtTopTagPlusTau23>=9"  ;;
    esac;

    case $HT in
    HT0) GO="${GO} -R ht500 ht5001000 htJet30j>500&&htJet30j<1000"  ;;
    HT1) GO="${GO} -R ht500 ht1000Inf htJet30j>=1000"  ;;
    HTDPhi) GO="${GO} -R ht500 ht1000Inf htJet30j>=1000 -R dp1 dp05 fabs(DeltaPhiLepW)>0.5 "  ;;
    HTStop) GO="${GO} -R ht500 ht1000Inf htJet30j>=1000 -R dp1 dp05 fabs(DeltaPhiLepW)>0.5 -A dp1 stopness (TopVarsMETovTopMin[0]-0.5)/0.5+(TopVarsMtopMin[0]-175)/175>1.25"  ;;
    HTTop) GO="${GO} -R ht500 ht1000Inf htJet30j>=1000 -R dp1 dp05 fabs(DeltaPhiLepW)>0.5 -A dp1 stopness (TopVarsMETovTopMin[0]-0.5)/0.5+(TopVarsMtopMin[0]-175)/175>1.25&&Topness>5"  ;;
    HTLowLepPt) GO="${GO} -R ht500 ht1000Inf htJet30j>=1000 -R 1tl 1tllowpt nTightLeps25==1&&LepGood_pt[tightLeps25idx[0]]<=25  -R dp1 dp00 fabs(DeltaPhiLepW)>0.0 -A dp1 stopness (TopVarsMETovTopMin[0]-0.5)/0.5+(TopVarsMtopMin[0]-175)/175>1.25&&Topness>5"  ;;
    HTLowLepPtDPhi) GO="${GO} -R ht500 ht1000Inf htJet30j>=1000 -R 1tl 1tllowpt nTightLeps25==1&&LepGood_pt[tightLeps25idx[0]]<=25"  ;;
    HTTTYes) GO="${GO} -R ht500 ht1000Inf htJet30j>=1000&&nHighPtTopTagPlusTau23>=1"  ;;
    HTTTNo) GO="${GO} -R ht500 ht1000Inf htJet30j>=1000&&nHighPtTopTagPlusTau23==0"  ;;
    esac;

    if [[ "$PRETEND" == "1" ]]; then
        echo "making datacard $OUT from makeShapeCardsSusy.py mca-Phys14_1l.txt bins/1l_CardsFullCutFlow.txt \"$EXPR\" \"$BINS\" $SYSTS $GO --dummyYieldsForZeroBkg;"
    else
        echo "making datacard $OUT from makeShapeCardsSusy.py mca-Phys14_1l.txt bins/1l_CardsFullCutFlow.txt \"$EXPR\" \"$BINS\" $SYSTS $GO --dummyYieldsForZeroBkg;"
        python makeShapeCardsSusy.py mca-Phys14_1l.txt bins/1l_CardsFullCutFlow.txt "$EXPR" "$BINS" $SYSTS -o $OUT $GO --dummyYieldsForZeroBkg;
        echo "  -- done at $(date)";
    fi;
}


function combineCardsSmart {
    DummyC=0
    AllC=0
    CMD=""
    for C in $*; do
        # missing datacards 
        test -f $C || continue;

	if grep -q "DummyContent" $C; then
	    echo "empty bin ${C}" >&2 
	    DummyC=$(($DummyC+1))
	    if grep -q "observation 0.0$" $C; then
		echo "this is not the way it was intended..."
	    fi
	fi

        grep -q "observation 0.0$" $C && continue
	AllC=$((AllC+1))
        CMD="${CMD} $(basename $C .card.txt)=$C ";
    done
    if [[ "$CMD" == "" ]]; then
        echo "Not any card found in $*" 1>&2 ;
    else
#	echo "combineCards.py $CMD" >&2
        combineCards.py $CMD
    fi
    if [[ DummyC>=0 ]]; then
	echo "In total $DummyC out of $AllC are empty, but taken into account by adding DummyContent." >&2
    fi
}

if [[ "$1" == "--pretend" ]]; then
    PRETEND=1; shift;
fi;
if [[ "$1" == "2lss-2012" ]]; then
    OPTIONS=" $OPTIONS -F sf/t $T/1_susyVars_2lssInc_v0/evVarFriend_{cname}.root "
    SYSTS="syst/susyDummy.txt"
    CnC_expr="1+4*(met_pt>120)+(htJet40j>400)+2*(nJet40>=4)"
    CnC_bins="[0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5]"
    MOD=inc;

    echo "Making individual datacards"
    for LL in ee em mm; do for LPt in 2020; do for SR in 0X 1X 2X+; do
        echo " --- CnC2012_${SR}_${LL} ---"
        #makeCard_2lss $CnC_expr $CnC_bins $SYSTS CnC2012_${SR}_${LL} "$OPTIONS";
    done; done; done
    echo "Making combined datacards"
    for D in $OUTDIR/T[0-9]*; do
        test -f $D/CnC2012_0X_ee.card.txt || continue
        (cd $D;
            for SR in 0X 1X 2X+; do
                combineCards.py CnC2012_${SR}_{ee,em,mm}.card.txt >  CnC2012_${SR}.card.txt
            done
            combineCards.py CnC2012_{0X,1X,2X+}.card.txt >  CnC2012.card.txt
        );
        echo "Made combined card $D/CnC2012.card.txt"
    done
    echo "Done at $(date)";

elif [[ "$1" == "2lss-2015" ]]; then
    OPTIONS=" $OPTIONS -F sf/t $T/1_susyVars_2lssInc_v2/evVarFriend_{cname}.root "
    SYSTS="syst/susyDummy.txt"
    CnC_expr="1+4*(met_pt>120)+(htJet40j>400)+2*(nJet40>=4)"
    CnC_bins="[0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5]"
    MOD=inc;

    echo "Making individual datacards"
    for LL in ee em mm; do for LPt in hh hl ll; do for SR in 0X 1X 2X 3X 2X+; do
    #for LL in ee em mm; do for LPt in hh hl ll ; do for SR in 0Xs 1Xs 2Xs 3Xs; do
        echo " --- CnC2015_${SR}_${LL}_${LPt} ---"
        makeCard_2lss $CnC_expr $CnC_bins $SYSTS CnC2015_${SR}_${LL}_${LPt} "$OPTIONS";
    done; done; done
    #exit
    echo "Making combined datacards"
    for D in $OUTDIR/T[0-9]*; do
        test -f $D/CnC2015_0X_ee_hh.card.txt || continue
        (cd $D && echo "    $D";
        for SR in 0X 1X 2X 3X 2X+; do
        #for SR in 0Xs 1Xs 2Xs 3Xs; do
            combineCardsSmart CnC2015_${SR}_{ee,em,mm}_hh.card.txt >  CnC2015_${SR}_hh.card.txt
            combineCardsSmart CnC2015_${SR}_{ee,em,mm}_{hh,hl,ll}.card.txt >  CnC2015_${SR}.card.txt
        done
        combineCardsSmart CnC2015_{0X,1X,2X+}.card.txt   >  CnC2015_2b.card.txt
        combineCardsSmart CnC2015_{0X,1X,2X+}_hh.card.txt   >  CnC2015_2b_hh.card.txt
        combineCardsSmart CnC2015_{0X,1X,2X,3X}_hh.card.txt >  CnC2015_3b_hh.card.txt
        combineCardsSmart CnC2015_{0X,1X,2X,3X}.card.txt >  CnC2015_3b.card.txt
        #combineCardsSmart CnC2015_{0Xs,1Xs,2Xs,3Xs}_hh.card.txt >  CnC2015_3bs_hh.card.txt
        #combineCardsSmart CnC2015_{0Xs,1Xs,2Xs,3Xs}.card.txt >  CnC2015_3bs.card.txt
        )
    done
    echo "Done at $(date)";

elif [[ "$1" == "2lss-2015x" ]]; then
    OPTIONS=" $OPTIONS -F sf/t $T/1_susyVars_2lssInc_v2/evVarFriend_{cname}.root "
    SYSTS="syst/susyDummy.txt"
    CnC_expr="1+4*(met_pt>120)+(htJet40j>400)+2*(nJet40>=4)"
    CnC_bins="[0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5]"
    MOD=excl;

    echo "Making individual datacards"
    for LL in ee em mm 3l; do for LPt in hh hl ll; do for SR in 0X 1X 2X 3X; do
        echo " --- CnC2015X_${SR}_${LL}_${LPt} ---"
        makeCard_2lss $CnC_expr $CnC_bins $SYSTS CnC2015X_${SR}_${LL}_${LPt} "$OPTIONS";
    done; done; done
    #exit
    echo "Making combined datacards"
    for D in $OUTDIR/T[0-9]*; do
        test -f $D/CnC2015X_0X_ee_hh.card.txt || continue
        (cd $D && echo "    $D";
        for SR in 0X 1X 2X 3X; do
            combineCardsSmart CnC2015X_${SR}_{ee,em,mm}_hh.card.txt >  CnC2015X_${SR}_hh.card.txt
            combineCardsSmart CnC2015X_${SR}_{ee,em,mm}_{hh,hl,ll}.card.txt >  CnC2015X_${SR}.card.txt
            combineCardsSmart CnC2015X_${SR}_{ee,em,mm,3l}_hh.card.txt >  CnC2015X_${SR}_hh_w3l.card.txt
            combineCardsSmart CnC2015X_${SR}_{ee,em,mm,3l}_{hh,hl,ll}.card.txt >  CnC2015X_${SR}_w3l.card.txt
        done
        combineCardsSmart CnC2015X_{0X,1X,2X,3X}_hh.card.txt >  CnC2015X_3b_hh.card.txt
        combineCardsSmart CnC2015X_{0X,1X,2X,3X}.card.txt >  CnC2015X_3b.card.txt
        combineCardsSmart CnC2015X_{0X,1X,2X,3X}_hh_w3l.card.txt >  CnC2015X_3b_hh_w3l.card.txt
        combineCardsSmart CnC2015X_{0X,1X,2X,3X}_w3l.card.txt >  CnC2015X_3b_w3l.card.txt
        )
    done
    echo "Done at $(date)";

fi

if [[ "$2" == "1l-2015" ]]; then
#    OPTIONS=" $OPTIONS -F sf/t $T/PHYS14_V2_Friends_MarkusHenningMerge_10GeVTightLepsInsteadOf25GeV/evVarFriend_{cname}.root "
    OPTIONS=" $OPTIONS -F sf/t $T/PHYS14_V3_Friends/evVarFriend_{cname}.root "
    SYSTS="syst/susyDummy.txt"
    CnC_expr="1" #not used as of now
    CnC_bins="[0.5,1.5]"


    echo "Making individual datacards"
    for ST in ST0 ST1 ST2 ST3 ST4; do for nJ in 45j 68j 6Infj 9Infj 68TTj 9InfTTj; do for nB in 0B 1B 2B 2Btop 3p; do for HT in HT0 HT1 HTDPhi HTStop HTTop HTLowLepPt HTLowLepPtDPhi HTTTYes HTTTNo; do
#    for ST in ST0; do for nJ in 45j 68j 6Infj 9Infj 68TTj 9InfTTj; do for nB in 2B 2Btop 3p; do for HT in HT0 HT1 HTDPhi HTStop HTTop HTLowLepPt HTTTYes HTTTNo; do
        echo " --- CnC2015X_${nB}_${ST}_${nJ}_${HT} ---"
        makeCard_1l $CnC_expr $CnC_bins $SYSTS CnC2015X_${nB}_${ST}_${nJ}_${HT} "$OPTIONS";
    done; done; done; done
    #exit
    echo "Making combined datacards"
    for D in $OUTDIR/T[0-9]*; do
        test -f $D/CnC2015X_2B_ST0_45j_HT0.card.txt || continue
        (cd $D && echo "    $D";
#        for nB in 0B 1B 2B 3p; do
        for nB in 2B 2Btop 3p; do
#	    echo "combineCardsSmart CnC2015X_${nB}_{ST0,ST1,ST2,ST3,ST4}_6Infj_{HT0,HT1}.card.txt >  CnC2015X_${nB}_standardnJ.card.txt"
            combineCardsSmart CnC2015X_${nB}_{ST0,ST1,ST2,ST3,ST4}_6Infj_{HT0,HT1}.card.txt >  CnC2015X_${nB}_standardnJ.card.txt
            combineCardsSmart CnC2015X_${nB}_{ST0,ST1,ST2,ST3,ST4}_{68j,9Infj}_{HT0,HT1}.card.txt >  CnC2015X_${nB}_finenJ.card.txt

            combineCardsSmart CnC2015X_${nB}_{ST0,ST1,ST2,ST3,ST4}_6Infj_HT1.card.txt >  CnC2015X_${nB}_standardnJ_HT1.card.txt # only high HT
            combineCardsSmart CnC2015X_${nB}_{ST0,ST1,ST2,ST3,ST4}_6Infj_HTLowLepPt.card.txt >  CnC2015X_${nB}_standardnJ_HTLowLepPt.card.txt # only high HT
            combineCardsSmart CnC2015X_${nB}_{ST0,ST1,ST2,ST3,ST4}_6Infj_HTLowLepPtDPhi.card.txt >  CnC2015X_${nB}_standardnJ_HTLowLepPtDPhi.card.txt # only high HT
            combineCardsSmart CnC2015X_${nB}_{ST0,ST1,ST2,ST3,ST4}_6Infj_HTDPhi.card.txt >  CnC2015X_${nB}_standardnJ_HTDPhi.card.txt # only high HT
            combineCardsSmart CnC2015X_${nB}_{ST0,ST1,ST2,ST3,ST4}_6Infj_HTStop.card.txt >  CnC2015X_${nB}_standardnJ_HTStop.card.txt # only high HT
            combineCardsSmart CnC2015X_${nB}_{ST0,ST1,ST2,ST3,ST4}_6Infj_HTTop.card.txt >  CnC2015X_${nB}_standardnJ_HTTop.card.txt # only high HT
            combineCardsSmart CnC2015X_${nB}_{ST0,ST1,ST2,ST3,ST4}_6Infj_{HTTTYes,HTTTNo}.card.txt >  CnC2015X_${nB}_standardnJ_HTTTYesNo.card.txt # only high HT
            combineCardsSmart CnC2015X_${nB}_{ST0,ST1,ST2,ST3,ST4}_6Infj_HTTTYes.card.txt >  CnC2015X_${nB}_standardnJ_HTTTYes.card.txt # only high HT
            combineCardsSmart CnC2015X_${nB}_{ST0,ST1,ST2,ST3,ST4}_6Infj_HTTTNo.card.txt >  CnC2015X_${nB}_standardnJ_HTTTNo.card.txt # only high HT

            combineCardsSmart CnC2015X_${nB}_{ST0,ST1,ST2,ST3,ST4}_{68j,9Infj}_HT1.card.txt >  CnC2015X_${nB}_finenJ_HT1.card.txt # only high HT
            combineCardsSmart CnC2015X_${nB}_{ST0,ST1,ST2,ST3,ST4}_{68j,9Infj}_HTLowLepPt.card.txt >  CnC2015X_${nB}_finenJ_HTLowLepPt.card.txt # only high HT
            combineCardsSmart CnC2015X_${nB}_{ST0,ST1,ST2,ST3,ST4}_{68j,9Infj}_HTLowLepPtDPhi.card.txt >  CnC2015X_${nB}_finenJ_HTLowLepPtDPhi.card.txt # only high HT
            combineCardsSmart CnC2015X_${nB}_{ST0,ST1,ST2,ST3,ST4}_{68j,9Infj}_HTDPhi.card.txt >  CnC2015X_${nB}_finenJ_HTDPhi.card.txt # only high HT
            combineCardsSmart CnC2015X_${nB}_{ST0,ST1,ST2,ST3,ST4}_{68j,9Infj}_HTStop.card.txt >  CnC2015X_${nB}_finenJ_HTStop.card.txt # only high HT
            combineCardsSmart CnC2015X_${nB}_{ST0,ST1,ST2,ST3,ST4}_{68j,9Infj}_HTTop.card.txt >  CnC2015X_${nB}_finenJ_HTTop.card.txt # only high HT
            combineCardsSmart CnC2015X_${nB}_{ST0,ST1,ST2,ST3,ST4}_{68j,9Infj}_{HTTTYes,HTTTNo}.card.txt >  CnC2015X_${nB}_finenJ_HTTTYesNo.card.txt # only high HT
            combineCardsSmart CnC2015X_${nB}_{ST0,ST1,ST2,ST3,ST4}_{68j,9Infj}_HTTTYes.card.txt >  CnC2015X_${nB}_finenJ_HTTTYes.card.txt # only high HT
            combineCardsSmart CnC2015X_${nB}_{ST0,ST1,ST2,ST3,ST4}_{68j,9Infj}_HTTTNo.card.txt >  CnC2015X_${nB}_finenJ_HTTTNo.card.txt # only high HT

            combineCardsSmart CnC2015X_${nB}_{ST0,ST1,ST2,ST3,ST4}_{68TTj,9InfTTj}_HT1.card.txt >  CnC2015X_${nB}_finenJ_HT1TT.card.txt # only high HT
            combineCardsSmart CnC2015X_${nB}_{ST0,ST1,ST2,ST3,ST4}_{68TTj,9InfTTj}_HTStop.card.txt >  CnC2015X_${nB}_finenJ_HTTopTT.card.txt # only high HT
            combineCardsSmart CnC2015X_${nB}_{ST0,ST1,ST2,ST3,ST4}_{68TTj,9InfTTj}_HTStop.card.txt >  CnC2015X_${nB}_MarkusProposal.card.txt
        done
        combineCardsSmart CnC2015X_{2B,3p}_standardnJ.card.txt >  CnC2015X_standardnJ.card.txt
        combineCardsSmart CnC2015X_{2B,3p}_finenJ.card.txt >  CnC2015X_finenJ.card.txt
        combineCardsSmart CnC2015X_{2B,3p}_standardnJ_HT1.card.txt >  CnC2015X_standardnJ_HT1.card.txt
        combineCardsSmart CnC2015X_{2B,3p}_standardnJ_HTLowLepPt.card.txt >  CnC2015X_standardnJ_HTLowLepPt.card.txt
        combineCardsSmart CnC2015X_{2B,3p}_standardnJ_HTLowLepPtDPhi.card.txt >  CnC2015X_standardnJ_HTLowLepPtDPhi.card.txt
        combineCardsSmart CnC2015X_{2B,3p}_standardnJ_HTDPhi.card.txt >  CnC2015X_standardnJ_HTDPhi.card.txt
        combineCardsSmart CnC2015X_{2B,3p}_standardnJ_HTStop.card.txt >  CnC2015X_standardnJ_HTStop.card.txt
        combineCardsSmart CnC2015X_{2B,3p}_standardnJ_HTTop.card.txt >  CnC2015X_standardnJ_HTTop.card.txt
        combineCardsSmart CnC2015X_{2B,3p}_standardnJ_HTTTYesNo.card.txt >  CnC2015X_standardnJ_HTTTYesNo.card.txt
        combineCardsSmart CnC2015X_{2B,3p}_standardnJ_HTTTYes.card.txt >  CnC2015X_standardnJ_HTTTYes.card.txt
        combineCardsSmart CnC2015X_{2B,3p}_standardnJ_HTTTNo.card.txt >  CnC2015X_standardnJ_HTTTNo.card.txt

        combineCardsSmart CnC2015X_{2B,3p}_finenJ_HT1.card.txt >  CnC2015X_finenJ_HT1.card.txt
        combineCardsSmart CnC2015X_{2B,3p}_finenJ_HTLowLepPt.card.txt >  CnC2015X_finenJ_HTLowLepPt.card.txt
        combineCardsSmart CnC2015X_{2B,3p}_finenJ_HTLowLepPtDPhi.card.txt >  CnC2015X_finenJ_HTLowLepPtDPhi.card.txt
        combineCardsSmart CnC2015X_{2B,3p}_finenJ_HTDPhi.card.txt >  CnC2015X_finenJ_HTDPhi.card.txt
        combineCardsSmart CnC2015X_{2B,3p}_finenJ_HTStop.card.txt >  CnC2015X_finenJ_HTStop.card.txt
        combineCardsSmart CnC2015X_{2B,3p}_finenJ_HTTop.card.txt >  CnC2015X_finenJ_HTTop.card.txt
        combineCardsSmart CnC2015X_{2B,3p}_finenJ_HTTTYesNo.card.txt >  CnC2015X_finenJ_HTTTYesNo.card.txt
        combineCardsSmart CnC2015X_{2B,3p}_finenJ_HTTTYes.card.txt >  CnC2015X_finenJ_HTTTYes.card.txt
        combineCardsSmart CnC2015X_{2B,3p}_finenJ_HTTTNo.card.txt >  CnC2015X_finenJ_HTTTNo.card.txt

        combineCardsSmart CnC2015X_{2B,3p}_finenJ_HT1TT.card.txt >  CnC2015X_finenJ_HT1TT.card.txt
        combineCardsSmart CnC2015X_{2B,3p}_finenJ_HTTopTT.card.txt >  CnC2015X_finenJ_HTTopTT.card.txt
        combineCardsSmart CnC2015X_{2Btop,3p}_MarkusProposal.card.txt >  CnC2015X_MarkusProposal.card.txt

        combineCardsSmart CnC2015X_standardnJ_{HT1,HTLowLepPt}.card.txt >  CnC2015X_standardnJ_HighLowLepPt.card.txt
        combineCardsSmart CnC2015X_finenJ_{HT1,HTLowLepPt}.card.txt >  CnC2015X_finenJ_HighLowLepPt.card.txt
        )
    done
    echo "Done at $(date)";

fi

