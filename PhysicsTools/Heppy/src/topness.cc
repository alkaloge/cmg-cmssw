/*
 *  topness.cpp
 *
 *  created on 01/15/13. by Francesco Costanza/ Dirk Kruecker
 *  DESY
 *  based on code from M.L.Graesser and J.Shelton
 @article{Graesser:2012qy,
      author         = "Graesser, Michael L. and Shelton, Jessie",
      title          = "{Hunting Asymmetric Stops}",
      journal        = "Phys.Rev.Lett.",
      volume         = "111",
      pages          = "121802",
      doi            = "10.1103/PhysRevLett.111.121802",
      year           = "2013",
      eprint         = "1212.4495",
      archivePrefix  = "arXiv",
      primaryClass   = "hep-ph",
      reportNumber   = "LA-UR-12-26897",
      SLACcitation   = "%%CITATION = ARXIV:1212.4495;%%",
 }
 *
 *  performs topness analysis on a single event
 *  using the Nelder-Mead method
 *
 *
 *  selection: -presumes event contains
 *
 *                at least two jets passing jet pre-selection requirements
 *
 *                with at least one of these jets b-tagged
 *
 *                one lepton passing event pre-selection
 *
 *
 *
 *****compiling *****
 *
 * make
 * needs ROOT
 *
 ********************
 *
 */

#include<iostream>
#include<string>
#include<vector>
#include<algorithm>
#include<cmath>
#include<stdexcept> // needed for exception runtime_error
#include<stdio.h>
#include<stdlib.h>
#include<valarray>
#include<map>

using namespace std;

// following file needed for pseduo-random number generator
#include "TRandom2.h"

// following file needed for minimization of topness
#include "PhysicsTools/Heppy/interface/asa047.h"

// TOPNESS header file
#include "PhysicsTools/Heppy/interface/topness.h"

// Example
int main(int argc, char *argv[]){
    int t0=clock();

    // initialize test event
    TLorentzVector mySamplelep;     mySamplelep.SetPxPyPzE(67.2,-45.2,-86.1,120.);

    TLorentzVector mySampleMET;     mySampleMET.SetPxPyPzE(-394,127.,0.,0.);

    TLorentzVector mySamplebjet;    mySamplebjet.SetPxPyPzE(143.7,-203.1,140.5,350.);
    TLorentzVector mySampleqcdjet1; mySampleqcdjet1.SetPxPyPzE(-95,47.,65.,170.);
    TLorentzVector mySampleqcdjet2; mySampleqcdjet2.SetPxPyPzE(109.,-58.,85.,200.);
    TLorentzVector mySampleqcdjet3; mySampleqcdjet3.SetPxPyPzE(169.,132.,-230.,467.);

    // initialize topness
    Topness::Topness top;

    std::vector<TLorentzVector> jets;
    std::vector<float> bDiscs;

    jets.push_back(mySampleqcdjet1); bDiscs.push_back(0);
    jets.push_back(mySampleqcdjet2); bDiscs.push_back(0);
    jets.push_back(mySampleqcdjet3); bDiscs.push_back(0);
    jets.push_back(mySamplebjet);    bDiscs.push_back(1);

    // get topness
    double lnT=log( top.GetTopness( jets, bDiscs, mySamplelep,  mySampleMET)); // ln T value of event

    cout << "The topness value on this event is " << lnT << endl;

    int t=clock()-t0;
    cout << "It took the program " << t << " clicks, or " << ((float) t)/CLOCKS_PER_SEC << " seconds." << endl;

    return 0;
}

Topness::Topness::Topness(){
    sigmat=15.;
    sigmaW=5.;
    sigmas=1000;

    // Default discriminator values (can be modified directly)

    //CSV 8 TeV
    //  bDisc_TightWP = 0.679;
    //  bDisc_LooseWP = 0.244;

    //CMVA 13 TeV
    //   bDisc_TightWP = 0.762;
    //   bDisc_LooseWP = 0.244; //not known, actually...

    //CSV v2 (CSV-IVF)
    bDisc_TightWP = 0.941;
    bDisc_LooseWP = 0.423;
    //bDisc_MediumWP = 0.814; // isn't used

    del=20.;
    lMAX=10;
}

void Topness::Topness::Sort( vector<const TLorentzVector*>& jets)
{
    sort( jets.begin(), jets.end(), ptSortPtr);
}

double Topness::Topness::GetTopness( const vector<TLorentzVector>& jets,
                                     const vector<float>& bDiscs,
                                     const TLorentzVector& lep,
                                     const TLorentzVector& MET)
{
    if ( jets.size() != bDiscs.size()){
        cout<<"Topness: jets.size() != BDiscs.size()"<<endl;
        return -1;
    }

    vector<const TLorentzVector*> noBJets;
    vector<const TLorentzVector*>  lBJets;
    vector<const TLorentzVector*>  tBJets;

    for ( unsigned int ijet = 0; ijet < jets.size(); ijet++){
        if (bDiscs.at(ijet) < bDisc_LooseWP) noBJets.push_back(&jets.at(ijet));
        else if (bDiscs.at(ijet) < bDisc_TightWP) lBJets.push_back(&jets.at(ijet));
        else tBJets.push_back(&jets.at(ijet));
    }

    Sort( noBJets);
    Sort(  lBJets);
    Sort(  tBJets);

    int nJets = jets.size();
    int nNoBJets = noBJets.size();
    int nLooseBJets = lBJets.size();
    int nTightBJets = tBJets.size();

    if ( nJets < 3)
        return -1.;

    vector<const TLorentzVector*> bjet;
    vector<const TLorentzVector*> jet;

    if ( nTightBJets == 0){
        if ( nLooseBJets == 0){
            bjet.push_back(noBJets.at(0));
            jet.push_back(noBJets.at(1));

            bjet.push_back(noBJets.at(0));
            jet.push_back(noBJets.at(2));

            bjet.push_back(noBJets.at(1));
            jet.push_back(noBJets.at(2));
        }
        else if ( nLooseBJets == 1){
            bjet.push_back(lBJets.at(0));
            jet.push_back(noBJets.at(0));

            bjet.push_back(lBJets.at(0));
            jet.push_back(noBJets.at(1));

            bjet.push_back(lBJets.at(0));
            jet.push_back(noBJets.at(1));
        }
        else{
            if ( nNoBJets > 0){
                bjet.push_back(lBJets.at(0));
                jet.push_back(lBJets.at(1));

                bjet.push_back(lBJets.at(0));
                jet.push_back(noBJets.at(0));

                bjet.push_back(lBJets.at(1));
                jet.push_back(noBJets.at(0));
            }
            else{
                bjet.push_back(lBJets.at(0));
                jet.push_back(lBJets.at(1));

                bjet.push_back(lBJets.at(0));
                jet.push_back(lBJets.at(2));

                bjet.push_back(lBJets.at(1));
                jet.push_back(lBJets.at(2));
            }
        }
    }
    else if ( nTightBJets == 1){
        if ( nLooseBJets == 0){
            bjet.push_back(tBJets.at(0));
            jet.push_back(noBJets.at(0));

            bjet.push_back(tBJets.at(0));
            jet.push_back(noBJets.at(1));

            bjet.push_back(tBJets.at(0));
            jet.push_back(noBJets.at(1));
        }
        else {
            if ( nNoBJets > 0){
                bjet.push_back(tBJets.at(0));
                jet.push_back(lBJets.at(0));

                bjet.push_back(tBJets.at(0));
                jet.push_back(noBJets.at(0));

                bjet.push_back(lBJets.at(0));
                jet.push_back(noBJets.at(0));
            }
            else{
                bjet.push_back(tBJets.at(0));
                jet.push_back(lBJets.at(0));

                bjet.push_back(tBJets.at(0));
                jet.push_back(lBJets.at(1));

                bjet.push_back(lBJets.at(0));
                jet.push_back(lBJets.at(1));
            }
        }
    }
    else{
        if ( nNoBJets > 0){
            bjet.push_back(tBJets.at(0));
            jet.push_back(tBJets.at(1));

            bjet.push_back(tBJets.at(0));
            jet.push_back(noBJets.at(0));

            bjet.push_back(tBJets.at(1));
            jet.push_back(noBJets.at(0));
        }
        else if ( nLooseBJets > 0){
            bjet.push_back(tBJets.at(0));
            jet.push_back(tBJets.at(1));

            bjet.push_back(tBJets.at(0));
            jet.push_back(lBJets.at(0));

            bjet.push_back(tBJets.at(1));
            jet.push_back(lBJets.at(0));
        }
        else{
            bjet.push_back(tBJets.at(0));
            jet.push_back(tBJets.at(1));

            bjet.push_back(tBJets.at(0));
            jet.push_back(tBJets.at(2));

            bjet.push_back(tBJets.at(1));
            jet.push_back(tBJets.at(2));
        }
    }

    /*if (nTightBJets == 0)
      return -1;

      TLorentzVector bjet, jet1, jet2;

      bjet = *tBJets.at(0);

      if ((nTightBJets>1)) {
      jet1=*tBJets.at(1);
      } else {
      jet1=*noBJets.at(0);
      if (nJets>1) {
      jet2=*noBJets.at(1);
      }
      }

      double topness_min = this->Calc( bjet, jet1, lep, MET);

      if (jet2.Pt() > 10.)
      {
      double topness_new = this->Calc( bjet, jet2, lep, MET);

      if (topness_new < topness_min)
      topness_min = topness_new;
      }

      return topness_min;
    */

    double topness = 9.e99;
    double topness_tmp = 9.e99;
    for ( unsigned int ijet = 0; ijet < jet.size(); ijet++){
        topness_tmp = this->Calc( *bjet.at(ijet), *jet.at(ijet), lep, MET);
        if ( topness_tmp < topness) topness = topness_tmp;
    }

    return topness;
}

double Topness::Topness::Calc( const TLorentzVector& bjet, const TLorentzVector& jet, const TLorentzVector& lep, const TLorentzVector& met){
    double pa[4];
    bjet.GetXYZT(pa);
    vector<double> b(pa, pa + 4 );
    jet.GetXYZT(pa);
    vector<double> j(pa, pa + 4 );
    lep.GetXYZT(pa);
    vector<double> l(pa, pa + 4 );
    met.GetXYZT(pa);
    vector<double> m(pa, pa + 4 );

    TopStatMET topstat( b, j, l, m, sigmat, sigmaW, sigmas, del, lMAX);

    double minpoint[4];

    return topstat.TopStatmin(minpoint);
};

bool Topness::ptSort(TLorentzVector p1, TLorentzVector p2){
    double pT1=p1.Pt();
    double pT2=p2.Pt();
    return (pT1>pT2);
}
bool Topness::ptSortPtr( const TLorentzVector* p1, const TLorentzVector* p2){
    double pT1=p1->Pt();
    double pT2=p2->Pt();
    return (pT1>pT2);
}

Topness::TopStatMET::TopStatMET(vector<double> ppb1, vector<double> ppb2, vector<double> ppl,vector<double> ppMET,
                                double sigmatt, double sigmaWW, double sigmass, double delp, int lMAXp)
    :pb1(ppb1), pb2(ppb2), pl(ppl), pMET(ppMET),
     sigmat(sigmatt), sigmaW(sigmaWW), sigmas(sigmass), del(delp), lMAX(lMAXp),
     g1(pb1,pb2,pl,pMET,sigmat,sigmaW,sigmas),
     g2(pb2,pb1,pl,pMET,sigmat,sigmaW,sigmas){};
// note the function gfuncMET is overloaded
// initialize g1 to physically measured momenta pb1, pb2, pl, pMET and to parameters of the statistic
// switch pb1 <-> pb2 to initialize g2

// this member function does all the work
// return min value
double Topness::TopStatMET::TopStatmin(double minpoint[])
{

    TRandom2 myranTopStatMET(7);

    // tolerance of Nelder-Mead method
    double ftol=0.000004;

    // initialize to some large number
    double finalmin=10000000000000.;

    // now make lMAX calls to Nelder-Mead
    for (int l=0;l<lMAX; l++)
    {

        double rndm=2*myranTopStatMET.Rndm()-1.;

        // starting point
        // attempt to choose starting point with adequeate coverage of phase space. 4000.=4 TeV.
        double startarr[]= {4000.*rndm,4000.*rndm,4000.*rndm,4000.*rndm};
        // defines start simplex
        double step[]={100,100,100,100};

        double fmin;

        // coordinates of minimum
        double min1[4],min2[4];
        int konvge = 10;
        int kcount = 500;
        int numres;
        int icount;
        int ifault;

        // will be value of function at minimum
        // first initialize at start point
        double fmin1 = g1(startarr);
        // call Nelder-Mead, using g1
        nelmin(g1,4,startarr,min1,&fmin1,ftol,step,konvge, kcount, &icount, &numres, &ifault );

        // now do the other combination
        double fmin2 = g2(startarr);
        nelmin(g2,4,startarr,min2,&fmin2,ftol,step,konvge, kcount, &icount, &numres, &ifault );

        // now pick smaller of two minimum
        if (fmin1<fmin2)
            fmin=fmin1;
        else
            fmin=fmin2;

        // if new minimum is smaller than current candidate global minimum
        // then it is the new candidate global minimum and store
        if (fmin < finalmin) {
            finalmin=fmin;
            if(fmin1<fmin2) {
                minpoint[0]=min1[0];
                minpoint[1]=min1[1];
                minpoint[2]=min1[2];
                minpoint[3]=min1[3];
            } else {
                minpoint[0]=min2[0];
                minpoint[1]=min2[1];
                minpoint[2]=min2[2];
                minpoint[3]=min2[3];
            }
        }

    }


    return finalmin;
};

double Topness::sign ( double x) {

    if (x>0)
        return  1.0;
    else    return -1.0;
};
