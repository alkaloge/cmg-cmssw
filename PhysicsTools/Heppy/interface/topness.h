#ifndef _TOPNESS_H_
#define _TOPNESS_H_

/*
 *  topness.h
 *
 *
 *  created 2/11/13.
 *
 */

#include<map>
#include<vector>
#include "TLorentzVector.h"

#include "asa047.h"

using namespace std;

namespace Topness
{
    class Topness
    {
    public:
        Topness();
        double GetTopness( const vector<TLorentzVector>& jets,
                           const vector<float>& bDiscs,
                           const TLorentzVector& lep,
                           const TLorentzVector& MET);

        double Calc( const TLorentzVector&,
                     const TLorentzVector&,
                     const TLorentzVector&,
                     const TLorentzVector&);

        double sigmat;
        double sigmaW;
        double sigmas;

        double bDisc_LooseWP;
        double bDisc_TightWP;

    private:
        void Sort( vector<const TLorentzVector*>& jets);
        vector<double> GetCoordinates( const TLorentzVector&);

        double del;       // step size in Nelder-Mead method (GeV)

        int lMAX;       // lMAX = number of random calls to Nelder-Mead method
    };

    bool ptSort(TLorentzVector  p1, TLorentzVector  p2);   // compares pT of p1 and p2
    bool ptSortPtr( const TLorentzVector* p1, const TLorentzVector* p2);   // compares pT of p1 and p2

    // BEGIN TOPNESS
    //
    // struct gfuncMET returns topness statistic
    // class TopStatMET performs the minimization of topness and calls gfuncMET

    // struct gfunctMET
    // initialized by b-jet, l, second jet  and MET momenta. Unknowns are 3-vec of neutrino and pZ component of the W.
    // using the measured transverse MET to solve for pWx, pWy. On-shell mass-conditions for v and W assumed.

    struct gfuncMET
    {
        unsigned int n;
        double mt,mW,mt2,mW2;
        //      static double mt,mW,mt2,mW2;

        // minkowski product
        inline double Lp(const vector<double> a, const vector<double> b){
            return a[3]*b[3]-a[0]*b[0]-a[1]*b[1]-a[2]*b[2];
        };
        // in place sum
        inline void ipS(vector<double>& a, const vector<double> b){
            a[0]+=b[0];
            a[1]+=b[1];
            a[2]+=b[2];
            a[3]+=b[3];
        };
        vector <double> pb2;
        vector <double> pb2l;
        vector <double> pb1, pb1b2l;
        vector <double> pl, pMET;


        //    vector<double> pb1, pb2, pl, pMET;
        //    vector<double> pb2l,pb1b2l;
        double sigmat4, sigmaW4, sigmas4;
        double mb12,ml2,mb2l2,mb1b2l2;
        double ma,mb,mc,md;
    gfuncMET(const vector<double> ppb1, const vector<double> ppb2, const vector<double> ppl, const vector<double> ppMET, double sigmat, double sigmaW, double sigmas)
    : pb2 (ppb2), pb2l(pb2), pb1 (ppb1), pb1b2l(pb1), pl (ppl), pMET (ppMET)
            {
                mt=172;     // top quark mass
                mW=80.4;    // W mass
                mt2=mt*mt;
                mW2=mW*mW;

                // initialize variables for quicker calculation
                sigmat4=pow(sigmat,4)/4;
                sigmaW4=pow(sigmaW,4)/4;
                sigmas4=pow(sigmas,4);
                mb12=Lp(pb1,pb1);
                ml2 =Lp(pl,pl);
                ipS(pb2l,pl);
                mb2l2=Lp(pb2l,pb2l);
                ipS(pb1b2l,pb2l);
                mb1b2l2=Lp(pb1b2l,pb1b2l);

                ma=(mb12+mW2-mt2)/2;
                mb=(mb2l2-mt2)/2;
                mc=(ml2-mW2)/2;
                md=4*mt2;

            }

        double operator()(double points[]) {  // points[0]=pv_x, points[1]=pv_y, points[2]=pv_z, points[3]=pW_z

            // pv_x, pv_y, pv_z
            double pvx=points[0];
            double pvy=points[1];
            double pvz=points[2];
            // neutrino energy assuming mass-shell condition
            double Ev=sqrt(pow(pvx,2)+pow(pvy,2)+pow(pvz,2));

            double pvarray[]={pvx,pvy,pvz,Ev};
            vector<double> pv(pvarray,pvarray+4);

            // pW_z
            double pWz=points[3];
            // W momenta from neutrino and MET
            double pWarray[]={-pvx+pMET[0],-pvy+pMET[1],pWz,
                              sqrt(pow(-pvx+pMET[0],2)+pow(-pvy+pMET[1],2)+pow(pWz,2)+mW2)};

            vector<double> pW(pWarray,pWarray+4);

            vector<double> pb1b2lWv(pWarray,pWarray+4);
            ipS(pb1b2lWv,pv);
            ipS(pb1b2lWv,pb1b2l);

            return pow(Lp(pb1,pW)+ma,2)/sigmat4
                +pow(Lp(pb2l,pv)+mb,2)/sigmat4
                +pow(Lp(pl,pv)+mc,2)/sigmaW4
                +pow(Lp(pb1b2lWv,pb1b2lWv)-md,2)/sigmas4;
        }

    };
    //double gfuncMET::mt=172;     // top quark mass
    //double gfuncMET::mW=80.4;    // W mass
    //double gfuncMET::mt2=gfuncMET::mt*gfuncMET::mt;
    //double gfuncMET::mW2=gfuncMET::mW*gfuncMET::mW;


    class TopStatMET
    {
        // perform minimization of the function gfuncMET using Nelder-Mead method
        // with lMAX calls of step-size del

    private:
        vector<double> pb1,pb2;
        vector <double> pl,pMET;
        double sigmat,sigmaW,sigmas,del;
        int lMAX;
        gfuncMET g1,g2;

    public:
        TopStatMET(vector<double>,vector<double>,vector<double>,vector<double>,double, double, double, double, int);
        //constructor
        double TopStatmin(double*);  // this member function does all the work
    };

    double sign ( double x);
}
#endif
