#ifndef FUNCTIONS_H
#define FUNCTIONS_H

// #include <stdio.h>
// #include <stdlib.h>
#include <iostream>
#include <cstdlib> //as stdlib.h                 
#include <cstdio>
#include <map>
#include <string>
#include <cmath>
#include "TH2F.h"
#include "TVector2.h"
#include "Math/GenVector/LorentzVector.h"
#include "Math/GenVector/PtEtaPhiM4D.h"

using namespace std;

//// UTILITY FUNCTIONS NOT IN TFORMULA ALREADY

float deltaPhi(float phi1, float phi2) {
    float result = phi1 - phi2;
    while (result > float(M_PI)) result -= float(2*M_PI);
    while (result <= -float(M_PI)) result += float(2*M_PI);
    return result;
}

float if3(bool cond, float iftrue, float iffalse) {
    return cond ? iftrue : iffalse;
}

float deltaR2(float eta1, float phi1, float eta2, float phi2) {
    float deta = std::abs(eta1-eta2);
    float dphi = deltaPhi(phi1,phi2);
    return deta*deta + dphi*dphi;
}
float deltaR(float eta1, float phi1, float eta2, float phi2) {
    return std::sqrt(deltaR2(eta1,phi1,eta2,phi2));
}

float Hypot(float x, float y) {
  return hypot(x,y);
}

float pt_2(float pt1, float phi1, float pt2, float phi2) {
    phi2 -= phi1;
    return hypot(pt1 + pt2 * std::cos(phi2), pt2*std::sin(phi2));
}

float mt_2(float pt1, float phi1, float pt2, float phi2) {
    return std::sqrt(2*pt1*pt2*(1-std::cos(phi1-phi2)));
}

float mass_2_ene(float ene1, float eta1, float phi1, float m1, float ene2, float eta2, float phi2, float m2) {
    typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > PtEtaPhiMVector;
    PtEtaPhiMVector unitp41(1.0,eta1,phi1,m1);
    PtEtaPhiMVector unitp42(1.0,eta2,phi2,m2);
    double theta1 = unitp41.Theta();
    double theta2 = unitp42.Theta();
    double pt1 = ene1*fabs(sin(theta1));
    double pt2 = ene2*fabs(sin(theta2));
    PtEtaPhiMVector p41(pt1,eta1,phi1,m1);
    PtEtaPhiMVector p42(pt2,eta2,phi2,m2);
    return (p41+p42).M();
}

float mass_2(float pt1, float eta1, float phi1, float m1, float pt2, float eta2, float phi2, float m2) {
    typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > PtEtaPhiMVector;
    PtEtaPhiMVector p41(pt1,eta1,phi1,m1);
    PtEtaPhiMVector p42(pt2,eta2,phi2,m2);
    return (p41+p42).M();
}

#include "TRandom3.h"
TRandom3 *randy = NULL;

float mass_2_smeared(float pt1, float eta1, float phi1, float m1, float pt2, float eta2, float phi2, float m2, int isData) {
    float finalpt1;
    float finalpt2;
    if (isData){
        finalpt1 = pt1;
        finalpt2 = pt2;
    }
    else {
        if (!randy) randy = new TRandom3(42);
        finalpt1 = pt1*(1.+0.34/pt1/1.4142*randy->Gaus(0.,1.) ) - 0.06957/pt1;
        finalpt2 = pt2*(1.+0.34/pt2/1.4142*randy->Gaus(0.,1.) ) - 0.06957/pt2;
    }

//    std::cout << "initial pT1: " << pt1 << " corrected pT1: " << finalpt1 << std::endl;
//    std::cout << "initial pT2: " << pt2 << " corrected pT2: " << finalpt2 << std::endl;

    typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > PtEtaPhiMVector;
    PtEtaPhiMVector p41(finalpt1,eta1,phi1,m1);
    PtEtaPhiMVector p42(finalpt2,eta2,phi2,m2);

    PtEtaPhiMVector p43(pt1,eta1,phi1,m1);
    PtEtaPhiMVector p44(pt2,eta2,phi2,m2);

    float finalmll = (p41+p42).M();
    float initialm = (p43+p44).M();

    //std::cout << "initial mll " << initialm << " final mll " << finalmll << std::endl;


    return (p41+p42).M();
}

float eta_2(float pt1, float eta1, float phi1, float m1, float pt2, float eta2, float phi2, float m2) {
    typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > PtEtaPhiMVector;
    PtEtaPhiMVector p41(pt1,eta1,phi1,m1);
    PtEtaPhiMVector p42(pt2,eta2,phi2,m2);
    return (p41+p42).Eta();
}

float pt_3(float pt1, float phi1, float pt2, float phi2, float pt3, float phi3) {
    phi2 -= phi1;
    phi3 -= phi1;
    return hypot(pt1 + pt2 * std::cos(phi2) + pt3 * std::cos(phi3), pt2*std::sin(phi2) + pt3*std::sin(phi3));
}

float mass_3(float pt1, float eta1, float phi1, float m1, float pt2, float eta2, float phi2, float m2, float pt3, float eta3, float phi3, float m3) {
    typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > PtEtaPhiMVector;
    PtEtaPhiMVector p41(pt1,eta1,phi1,m1);
    PtEtaPhiMVector p42(pt2,eta2,phi2,m2);
    PtEtaPhiMVector p43(pt3,eta3,phi3,m3);
    return (p41+p42+p43).M();
}

float pt_4(float pt1, float phi1, float pt2, float phi2, float pt3, float phi3, float pt4, float phi4) {
    phi2 -= phi1;
    phi3 -= phi1;
    phi4 -= phi1;
    return hypot(pt1 + pt2 * std::cos(phi2) + pt3 * std::cos(phi3) + pt4 * std::cos(phi4), pt2*std::sin(phi2) + pt3*std::sin(phi3) + pt4*std::sin(phi4));
}
 
float mass_4(float pt1, float eta1, float phi1, float m1, float pt2, float eta2, float phi2, float m2, float pt3, float eta3, float phi3, float m3, float pt4, float eta4, float phi4, float m4) {
    typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > PtEtaPhiMVector;
    PtEtaPhiMVector p41(pt1,eta1,phi1,m1);
    PtEtaPhiMVector p42(pt2,eta2,phi2,m2);
    PtEtaPhiMVector p43(pt3,eta3,phi3,m3);
    PtEtaPhiMVector p44(pt4,eta4,phi4,m4);
    return (p41+p42+p43+p44).M();
}

float mt_llv(float ptl1, float phil1, float ptl2, float phil2, float ptv, float phiv) {
    float px = ptl1*std::cos(phil1) + ptl2*std::cos(phil2) + ptv*std::cos(phiv);
    float py = ptl1*std::sin(phil1) + ptl2*std::sin(phil2) + ptv*std::sin(phiv);
    float ht = ptl1+ptl2+ptv;
    return std::sqrt(std::max(0.f, ht*ht - px*px - py*py));
}

float mt_lllv(float ptl1, float phil1, float ptl2, float phil2, float ptl3, float phil3, float ptv, float phiv) {
    float px = ptl1*std::cos(phil1) + ptl2*std::cos(phil2) + ptl3*std::cos(phil3) + ptv*std::cos(phiv);
    float py = ptl1*std::sin(phil1) + ptl2*std::sin(phil2) + ptl3*std::sin(phil3) + ptv*std::sin(phiv);
    float ht = ptl1+ptl2+ptl3+ptv;
    return std::sqrt(std::max(0.f, ht*ht - px*px - py*py));
}


float mtw_wz3l(float pt1, float eta1, float phi1, float m1, float pt2, float eta2, float phi2, float m2, float pt3, float eta3, float phi3, float m3, float mZ1, float met, float metphi) 
{
    if (abs(mZ1 - mass_2(pt1,eta1,phi1,m1,pt2,eta2,phi2,m2)) < 0.01) return mt_2(pt3,phi3,met,metphi);
    if (abs(mZ1 - mass_2(pt1,eta1,phi1,m1,pt3,eta3,phi3,m3)) < 0.01) return mt_2(pt2,phi2,met,metphi);
    if (abs(mZ1 - mass_2(pt2,eta2,phi2,m2,pt3,eta3,phi3,m3)) < 0.01) return mt_2(pt1,phi1,met,metphi);
    return 0;
}

float mt_lu_cart(float lep_pt, float lep_phi, float u_x, float u_y)
{
    float lep_px = lep_pt*std::cos(lep_phi), lep_py = lep_pt*std::sin(lep_phi);
    float u = hypot(u_x,u_y);
    float uDotLep = u_x*lep_px + u_y*lep_py;
    return sqrt(2*lep_pt*sqrt(u*u+lep_pt*lep_pt+2*uDotLep) + 2*uDotLep + 2*lep_pt*lep_pt);
}

float u1_2(float met_pt, float met_phi, float ref_pt, float ref_phi) 
{
    float met_px = met_pt*std::cos(met_phi), met_py = met_pt*std::sin(met_phi);
    float ref_px = ref_pt*std::cos(ref_phi), ref_py = ref_pt*std::sin(ref_phi);
    float ux = - met_px + ref_px, uy = - met_py + ref_py;
    return (ux*ref_px + uy*ref_py)/ref_pt;
}
float u2_2(float met_pt, float met_phi, float ref_pt, float ref_phi)
{
    float met_px = met_pt*std::cos(met_phi), met_py = met_pt*std::sin(met_phi);
    float ref_px = ref_pt*std::cos(ref_phi), ref_py = ref_pt*std::sin(ref_phi);
    float ux = - met_px + ref_px, uy = - met_py + ref_py;
    return (ux*ref_py - uy*ref_px)/ref_pt;
}

float met_cal(float met_pt, float met_phi, float lep_pt, float lep_phi, float u_coeff, float u_syst)
{
    float met_px = met_pt*std::cos(met_phi), met_py = met_pt*std::sin(met_phi);
    float lep_px = lep_pt*std::cos(lep_phi), lep_py = lep_pt*std::sin(lep_phi);
    float ux = met_px + lep_px, uy = met_py + lep_py;
    float metcal_px = - u_coeff*ux*(1+u_syst) - lep_px, metcal_py = - u_coeff*uy*(1+u_syst) - lep_py;
    return hypot(metcal_px,metcal_py);
}

float _puw2016_nTrueInt_BF[60] = {0.0004627598152210959, 0.014334910915287028, 0.01754727657726197, 0.03181477917631854, 0.046128282569231016, 0.03929080994013006, 0.057066019809589925, 0.19570744862221007, 0.3720256062526554, 0.6440076202772811, 0.9218024454406528, 1.246743510634073, 1.5292543296414058, 1.6670061646418215, 1.7390553377117133, 1.6114721876895595, 1.4177294439817985, 1.420132866045718, 1.3157656415540477, 1.3365188060918483, 1.1191478126677334, 0.9731079434848392, 0.9219564145009487, 0.8811793391804676, 0.7627315352977334, 0.7265186492688713, 0.558602385324645, 0.4805954159733825, 0.34125298049234554, 0.2584848657646724, 0.1819638766151892, 0.12529545619337035, 0.11065705912071645, 0.08587356267495487, 0.09146322371620583, 0.11885517671051576, 0.1952483711863489, 0.23589115679998116, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
float puw2016_nTrueInt_BF(int nTrueInt) { if (nTrueInt<60) return _puw2016_nTrueInt_BF[nTrueInt]; else return 0; }

float _puw2016_nTrueInt_36fb[100] = {0.3505407355600995, 0.8996968628890968, 1.100322319466069, 0.9562526765089195, 1.0366251229154624, 1.0713954619016586, 0.7593488199769544, 0.47490309461978414, 0.7059895997695581, 0.8447022252423783, 0.9169159386164522, 1.0248924033173097, 1.0848877947714115, 1.1350984224561655, 1.1589888429954602, 1.169048420382294, 1.1650383018054549, 1.1507200023444994, 1.1152571438041776, 1.0739529436969637, 1.0458014000030829, 1.032500407707141, 1.0391236062781293, 1.041283620738903, 1.0412963370894526, 1.0558823002770783, 1.073481674823461, 1.0887053272606795, 1.1041701696801014, 1.123218903738397, 1.1157169321377927, 1.1052520327174429, 1.0697489590429388, 1.0144652740600584, 0.9402657069968621, 0.857142825520793, 0.7527112615290031, 0.6420618248685722, 0.5324755829715156, 0.4306470627563325, 0.33289171600176093, 0.24686361729094983, 0.17781595237914027, 0.12404411884835284, 0.08487088505600057, 0.056447805688061216, 0.03540829360547507, 0.022412461576677457, 0.013970541270658443, 0.008587896629717911, 0.004986410514292661, 0.00305102303701641, 0.001832072556146534, 0.0011570757619737708, 0.0008992999249003301, 0.0008241241729452477, 0.0008825716073180279, 0.001187003960081393, 0.0016454104270429153, 0.0022514113879764414, 0.003683196037880878, 0.005456695951503178, 0.006165248770884191, 0.007552675218762607, 0.008525338219226993, 0.008654690499815343, 0.006289068906974821, 0.00652551838513972, 0.005139581024893171, 0.005115751962934923, 0.004182527768384693, 0.004317593022028565, 0.0035749335962533355, 0.003773660372937113, 0.002618732319396435, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
float puw2016_nTrueInt_36fb(int nTrueInt) { if (nTrueInt<100) return _puw2016_nTrueInt_36fb[nTrueInt]; else return 0; }


float _puw_marc[100] = {0.436283166038, 0.918116564281, 1.33519020667, 0.951489790052, 1.08238493627, 1.21717271857, 0.79862426782, 0.48746564303, 0.744027463114, 0.885321584404, 0.974389321007, 1.06664838577, 1.12630370013, 1.17540434521, 1.20867880114, 1.20563961433, 1.2023653373, 1.181665256, 1.1452017341, 1.09507135359, 1.0691734481, 1.05074832463, 1.05194701757, 1.05146058305, 1.05001362966, 1.05339366189, 1.07098159855, 1.08754654701, 1.09596673184, 1.10279812058, 1.08881427644, 1.08075170669, 1.04354652211, 0.983387116092, 0.913947413337, 0.824243524918, 0.717414514035, 0.607106245317, 0.500460658887, 0.402180628848, 0.308931660162, 0.228191744417, 0.163558426749, 0.113614405206, 0.0775267958816, 0.0509487615982, 0.03193482641, 0.0202508940305, 0.0122188696444, 0.00748567745608, 0.00441765722763, 0.00264549673999, 0.00157754766081, 0.000965447543018, 0.000738368527738, 0.000659309925226, 0.000715517272062, 0.00094209777685, 0.00137951048957, 0.0018788996289, 0.00297177633742, 0.00401802341487, 0.00507326229236, 0.00546027560291, 0.00430719634822, 0.00596840555965, 0.00452179215829, 0.00480779513226, 0.00348281397548, 0.00338582312161, 0.00282389625724, 0.00254018692427, 0.001996654177, 0.00168910262341, 0.00169353843379, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
float puw_marc(int nTrueInt) { if (nTrueInt<100) return _puw_marc[nTrueInt]; else return 0.; }


// functions to assess if events pass given ID cuts
// isEB can be defined as (LepGood1_etaSc)<1.479 
// note that 2016 cut-based ID defines thesholds for EB and EE using SuperCluster eta
// the real ID WP part is in LepGood1_tightId, LepGood1_lostHits and LepGood1_convVeto
// dxy and dz are not part of the official ID WP, but we use the suggested thresholds anyway
// https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2
// 
// list of functions to manage IDs
// those marked with ** are work in progress (problems with string, format is n tcompatible with TTree::Draw() used by mcPlots.py)
// -------------------------------
// pass_dxy_dz
// pass_lostHits_conVe
// pass_looseIDnoIso_2016
// pass_mediumIDnoIso_2016
// pass_tightIDnoIso_2016
// pass_workingPointIDnoIso_2016 **
// pass_isolation_2016 **
// passFakerateNumerator2016 **
// isInFakerateApplicationRegion2016 **
// pass_isolation_WP
// pass_FakerateNumerator2016
// pass_FakerateApplicationRegion2016
//
//
// -------------------------------


// pass dxy and dz
bool pass_dxy_dz(const bool isEB = true, 
		 const float LepGood1_dxy = -999, 
		 const float LepGood1_dz = -999
		 ) 
{
  
  if (isEB) return (abs(LepGood1_dxy) < 0.05 && abs(LepGood1_dz) < 0.1);
  else      return (abs(LepGood1_dxy) < 0.1  && abs(LepGood1_dz) < 0.2);

}

// missing hits and conversion veto
bool pass_lostHits_conVeto(const int LepGood1_lostHits = -999, 
			   const int LepGood1_convVeto = -999
			   ) 
{
  return (LepGood1_lostHits <= 1 && LepGood1_convVeto == 1);
}


// loose ID no isolation
bool pass_looseIDnoIso_2016(const bool  isEB = true, 
			    const int   LepGood1_tightId = -1, 
			    const float LepGood1_dxy = -999, 
			    const float LepGood1_dz = -999,
			    const int   LepGood1_lostHits = -1,
			    const int   LepGood1_convVeto = -999
			    ) 
{

  return (LepGood1_tightId >= 1 && pass_dxy_dz(isEB,LepGood1_dxy,LepGood1_dz) && pass_lostHits_conVeto(LepGood1_lostHits,LepGood1_convVeto) );

}

// medium ID no isolation
bool pass_mediumIDnoIso_2016(const bool  isEB = true, 
			     const int   LepGood1_tightId = -1, 
			     const float LepGood1_dxy = -999, 
			     const float LepGood1_dz = -999,
			     const int   LepGood1_lostHits = -1,
			     const int   LepGood1_convVeto = -999
			     ) 
{

  return (LepGood1_tightId >= 2 && pass_dxy_dz(isEB,LepGood1_dxy,LepGood1_dz) && pass_lostHits_conVeto(LepGood1_lostHits,LepGood1_convVeto) );

}

// tight ID no isolation
bool pass_tightIDnoIso_2016(const bool  isEB = true, 
			    const int   LepGood1_tightId = -1, 
			    const float LepGood1_dxy = -999, 
			    const float LepGood1_dz = -999,
			    const int   LepGood1_lostHits = -1,
			    const int   LepGood1_convVeto = -999
			    ) 
{

  return (LepGood1_tightId >= 3 && pass_dxy_dz(isEB,LepGood1_dxy,LepGood1_dz) && pass_lostHits_conVeto(LepGood1_lostHits,LepGood1_convVeto) );

}

/////////////////////////////////////////////////////
//
// Following commented functions are work in progress
// Problems in using string
//
/////////////////////////////////////////////////////

// bool pass_workingPointIDnoIso_2016(const string workingPoint = "loose", // loose, medium, tight
// 				   const bool  isEB = true, 
// 				   const int   LepGood1_tightId = -1, 
// 				   const float LepGood1_dxy = -999, 
// 				   const float LepGood1_dz = -999,
// 				   const int   LepGood1_lostHits = -1,
// 				   const int   LepGood1_convVeto = -999
// 				   ) 
// {

//   if      (workingPoint == "loose" ) return pass_looseIDnoIso_2016(  isEB,LepGood1_tightId,LepGood1_dxy,LepGood1_dz,LepGood1_lostHits,LepGood1_convVeto);
//   else if (workingPoint == "medium") return pass_mediumIDnoIso_2016( isEB,LepGood1_tightId,LepGood1_dxy,LepGood1_dz,LepGood1_lostHits,LepGood1_convVeto);
//   else if (workingPoint == "tight" ) return pass_tightIDnoIso_2016(  isEB,LepGood1_tightId,LepGood1_dxy,LepGood1_dz,LepGood1_lostHits,LepGood1_convVeto);
//   else {
//     cout << "Error in function pass_workingPointIDnoIso_2016(): undefined working point "<< workingPoint << ", please check. Exiting ..." <<endl;
//     exit(EXIT_FAILURE);
//   }

// }


// bool pass_isolation_2016(const string workingPoint = "loose", // loose, medium, tight, custom
// 			 const bool   isEB = true,
// 			 const float  LepGood1_relIso04EA = -1,
// 			 )
// {

//   // WARNING: test that strings are accepted by mc*.py, currently they are not
//   // function format should be compatible with TTre::Draw()
//   std::map<string,float> workingPointIsolation_EB;
//   workingPointIsolation_EB["veto"  ] = 0.175; 
//   workingPointIsolation_EB["loose" ] = 0.0994; 
//   workingPointIsolation_EB["medium"] = 0.0695; 
//   workingPointIsolation_EB["tight" ] = 0.0588; 
//   workingPointIsolation_EB["custom"] = 0.2; 
//   std::map<std::string,float> workingPointIsolation_EE;
//   workingPointIsolation_EE["veto"  ] = 0.159; 
//   workingPointIsolation_EE["loose" ] = 0.107; 
//   workingPointIsolation_EE["medium"] = 0.0821; 
//   workingPointIsolation_EE["tight" ] = 0.0571; 
//   workingPointIsolation_EE["custom"] = 0.0821;

//   if (isEB) return LepGood1_relIso04EA < workingPointIsolation_EB[workingPoint];
//   else      return LepGood1_relIso04EA < workingPointIsolation_EE[workingPoint];

// }

// bool passFakerateNumerator2016(const string workingPoint = "loose", // loose, medium, tight
// 			       const bool   isEB = true, 
// 			       const int    LepGood1_tightId = -1, 
// 			       const float  LepGood1_dxy = -999, 
// 			       const float  LepGood1_dz = -999,
// 			       const int    LepGood1_lostHits = -1,
// 			       const int    LepGood1_convVeto = -999,
// 			       const float  LepGood1_relIso04EA = -1,
// 			       const bool   useCustomRelIso04EA = true // use user defined isolation threshold, not the E/gamma value
// 			       ) 
// {

//   return (pass_workingPointIDnoIso_2016(workingPoint,isEB,LepGood1_tightId,LepGood1_dxy,LepGood1_dz,LepGood1_lostHits,LepGood1_convVeto) 
// 	  && 
// 	  pass_isolation_2016(workingPoint,isEB,LepGood1_relIso04EA,useCustomRelIso04EA)
// 	  );

// }


// bool isInFakerateApplicationRegion2016(const string workingPoint = "loose", // loose, medium, tight
// 				       const bool   isEB = true, 
// 				       const int    LepGood1_tightId = -1, 
// 				       const float  LepGood1_dxy = -999, 
// 				       const float  LepGood1_dz = -999,
// 				       const int    LepGood1_lostHits = -1,
// 				       const int    LepGood1_convVeto = -999,
// 				       const float  LepGood1_relIso04EA = -1,
// 				       const bool   useCustomRelIso04EA = true // use user defined isolation threshold, not the E/gamma value
// 				       ) 
// {

//   return (not passFakerateNumerator2016(workingPoint,isEB,
// 					LepGood1_tightId,LepGood1_dxy,LepGood1_dz,LepGood1_lostHits,LepGood1_convVeto,
// 					LepGood1_relIso04EA,useCustomRelIso04EA
// 					)
// 	  );

// }

//==========================

bool pass_isolation_WP(const bool isEB = true, const float  LepGood1_relIso04EA = -1)
{
  // custom WP for 2016 data (before final Legacy ReReco, it could be changed for the new last one)
  return (LepGood1_relIso04EA < (isEB ? 0.2 : 0.0821)); // custom value for EB, medium WP for EE
}

//==========================

bool pass_looseIsolation_2016(const bool   isEB = true,
			      const float  LepGood1_relIso04EA = -1
			      )
{

  if (isEB) return LepGood1_relIso04EA < 0.0994;
  else      return LepGood1_relIso04EA < 0.107;

}

//==========================


bool pass_mediumIsolation_2016(const bool   isEB = true,
			       const float  LepGood1_relIso04EA = -1
			       )
{

  if (isEB) return LepGood1_relIso04EA < 0.0695;
  else      return LepGood1_relIso04EA < 0.0821;

}

//==========================


bool pass_tightIsolation_2016(const bool   isEB = true,
			      const float  LepGood1_relIso04EA = -1
			      )
{

  if (isEB) return LepGood1_relIso04EA < 0.0588;
  else      return LepGood1_relIso04EA < 0.0571;

}

//==========================

bool pass_FakerateNumerator_loose2016(const bool   isEB = true, 
				      const int    LepGood1_tightId = -1, 
				      const float  LepGood1_dxy = -999, 
				      const float  LepGood1_dz = -999,
				      const int    LepGood1_lostHits = -1,
				      const int    LepGood1_convVeto = -999,
				      const float  LepGood1_relIso04EA = -1
				      ) 
{
  
    return (pass_looseIDnoIso_2016(isEB,LepGood1_tightId,LepGood1_dxy,LepGood1_dz,LepGood1_lostHits,LepGood1_convVeto)
	    && 
	    pass_looseIsolation_2016(isEB,LepGood1_relIso04EA)
	    );

}

//============================================

bool pass_FakerateNumerator_medium2016(const bool   isEB = true, 
				       const int    LepGood1_tightId = -1, 
				       const float  LepGood1_dxy = -999, 
				       const float  LepGood1_dz = -999,
				       const int    LepGood1_lostHits = -1,
				       const int    LepGood1_convVeto = -999,
				       const float  LepGood1_relIso04EA = -1
				       ) 
{
  
    return (pass_mediumIDnoIso_2016(isEB,LepGood1_tightId,LepGood1_dxy,LepGood1_dz,LepGood1_lostHits,LepGood1_convVeto)
	    && 
	    pass_mediumIsolation_2016(isEB,LepGood1_relIso04EA)
	    );

}

//============================================

bool pass_FakerateNumerator2016(const bool   isEB = true, 
				const int    LepGood1_tightId = -1, 
				const float  LepGood1_dxy = -999, 
				const float  LepGood1_dz = -999,
				const int    LepGood1_lostHits = -1,
				const int    LepGood1_convVeto = -999,
				const float  LepGood1_relIso04EA = -1
				) 
{

  // EB, loose ID + iso < 0.2
  // EE full medium ID + iso

  if (isEB) {
    return (pass_looseIDnoIso_2016(isEB,LepGood1_tightId,LepGood1_dxy,LepGood1_dz,LepGood1_lostHits,LepGood1_convVeto)
	    && 
	    pass_isolation_WP(isEB,LepGood1_relIso04EA)
	    );
  } else {
    return (pass_mediumIDnoIso_2016(isEB,LepGood1_tightId,LepGood1_dxy,LepGood1_dz,LepGood1_lostHits,LepGood1_convVeto)
	    &&
	    pass_isolation_WP(isEB,LepGood1_relIso04EA)
	    );
  }

}

//==========================


bool pass_FakerateApplicationRegion2016(const bool   isEB = true, 
					const int    LepGood1_tightId = -1, 
					const float  LepGood1_dxy = -999, 
					const float  LepGood1_dz = -999,
					const int    LepGood1_lostHits = -1,
					const int    LepGood1_convVeto = -999,
					const float  LepGood1_relIso04EA = -1
					) 
{

  return (not pass_FakerateNumerator2016(isEB,LepGood1_tightId,LepGood1_dxy,LepGood1_dz,LepGood1_lostHits,LepGood1_convVeto,LepGood1_relIso04EA));

}


//==================================================

bool pass_FakerateNum_debug(const bool  isEB = true, 
			    const int   LepGood1_tightId = -1, 
			    const float LepGood1_dxy = -999, 
			    const float LepGood1_dz = -999,
			    const int   LepGood1_lostHits = -1,
			    const int   LepGood1_convVeto = -999,
			    const float LepGood1_relIso04EA = -1
			    ) 
{

  // EB, loose ID + iso < 0.2
  // EE full medium ID + iso

  if (isEB) {
    return (LepGood1_tightId >= 1 && abs(LepGood1_dxy) <= 0.05 && abs(LepGood1_dxy) <= 0.1 && LepGood1_lostHits <= 1 && LepGood1_convVeto == 1 && LepGood1_relIso04EA <= 0.2);
  } else {
    return (LepGood1_tightId >= 2 && abs(LepGood1_dxy) <= 0.1 && abs(LepGood1_dxy) <= 0.2 && LepGood1_lostHits <= 1 && LepGood1_convVeto == 1 && LepGood1_relIso04EA <= 0.0821);
  }

}
//==================================================

TVector2 tkmetEleCorr(float tkmet_pt, float tkmet_phi, float lep_pt, float lep_phi) {

  TVector2 trkmet_corr;
  trkmet_corr.SetMagPhi(tkmet_pt, tkmet_phi);

  // when the electron is not compatible with the primary vertex, its track is not used to compute tkMet (it is a bug in our ntuples)
  // in that case, add the electron back
  // We have the following (assuming vectorial object in the equation)
  // TkMEt_corr = -Sum(pT_tracks_noBadEle) - pT_badEle
  // in the ntuples we have TkMEt = -Sum(pT_tracks_noBadEle)
  
  // here we correct the tkMet 
  TVector2 badEle;      
  badEle.SetMagPhi(lep_pt,lep_phi);
  trkmet_corr -= badEle;  

  return trkmet_corr;

}

//==================================================    

float tkmetEleCorr_pt(float tkmet_pt, float tkmet_phi, float lep_pt, float lep_phi, bool eleTrackIsVertexCompatible) {

  if (eleTrackIsVertexCompatible) return tkmet_pt;

  TVector2 trkmet_corr = tkmetEleCorr(tkmet_pt, tkmet_phi, lep_pt, lep_phi);
  return trkmet_corr.Mod();

}


//==================================================


float tkmetEleCorr_phi(float tkmet_pt, float tkmet_phi, float lep_pt, float lep_phi, bool eleTrackIsVertexCompatible) {

  if (eleTrackIsVertexCompatible) return tkmet_phi;

  TVector2 trkmet_corr = tkmetEleCorr(tkmet_pt, tkmet_phi, lep_pt, lep_phi);
  return trkmet_corr.Phi();

}


//==================================================

// call like 
// tkmt_tkmetEleCorr(met_trkPt,
// 		     met_trkPhi,
// 		     ptElFull(LepGood1_calPt,LepGood1_eta),
// 		     LepGood1_phi, 
// 		     pass_dxy_dz(abs(LepGood1_eta)<1.479, LepGood1_dxy, LepGood1_dz) && pass_lostHits_conVeto(LepGood1_lostHits, LepGood1_convVeto))
// The correct definition is actually requiring |dz| < 0.1, which is the condition to be fullfilled to have a charge particle used for TrkMET in out ntuples

float tkmt_tkmetEleCorr(float tkmet_pt, float tkmet_phi, float lep_pt, float lep_phi, bool eleTrackIsVertexCompatible) {

  if (eleTrackIsVertexCompatible) {

    return mt_2(tkmet_pt, tkmet_phi, lep_pt, lep_phi);

  } else {

    // when the electron is not compatible with the primary vertex, its track is not used to compute tkMet (it is a bug in our ntuples)
    // in that case, add the electron back
    // We have the following (assuming vectorial object in the equation)
    // TkMEt_corr = -Sum(pT_tracks_noBadEle) - pT_badEle
    // in the ntuples we have TkMEt = -Sum(pT_tracks_noBadEle)

    TVector2 trkmet_corr = tkmetEleCorr(tkmet_pt, tkmet_phi, lep_pt, lep_phi); 
    return mt_2(trkmet_corr.Mod(),trkmet_corr.Phi(),lep_pt,lep_phi);

  }

}

//==================================================



void functions() {}

#endif
