/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   globals.hpp
 * Author: Hiromi WL Koh
 *
 * Last modified on 13 June, 2018, 4:00 PM
 */

#ifndef GLOBALS_HPP
#define GLOBALS_HPP

#include "geneClass.hpp"
#include "subjectClass.hpp"
#include "PAMClass.hpp"

#include <stdio.h>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

const double NA_VAL = 1e6;
const int NA_INT = 9999999;

struct Input_t{
   string dataZ,dataY,dataX;
   string withinNet,betweenNet;
   string pathwayfile,subtypefile, priorfile;
   string directory;
   string normalize_Zby = "Y";
   vector<int> cid_Z, cid_Y, cid_X;
   bool log_Z = true;
   bool log_Y= true;
   bool log_X = true;
   bool ztrans_Z = true;
   bool ztrans_Y = true;
   bool ztrans_X = true;

   bool analyzeZ = false;
   bool analyzeY = false;
   bool knnimpute = false;
   bool crossV = true;
   bool enrichment = true;
   bool prior_input = false;

   int blocksize =1500;
   int minobs=1; 
   int knnK=10; 
   int numFold = 10;
   double minprop=0.5;
   double minthres=NA_VAL;

   int minbgSize = 1;
   int minsigSize = 3;
   double bgProp = 0.5;

};


struct PAM_t{
    string geneA, geneB, identifier;
    string dataT;
    int NumFD_btw;
    int NumFD_within;
     vector<set<string> > FDneighbor;
    vector<double> prop_agree;
    vector<double> FDbtw_dik_bar, FDwithin_dik_bar;
    vector<double> d_ik, xbar_ik, mk;
    vector<double> dikstar_new;
    vector<string> classLab;
    double s_o,s_i, xbar_i;
    vector<int> GRPsize;
};



// General functions //
vector<string> stringSplit(string str, char delim);
vector<string> stringSplit2(string str, string delim);
string concatenate(vector<string> vec, char delim);
bool GREP(string STRING, string PATTERN);
vector<int> whichMin(vector<double> vec);
int whichMax(vector<double> vec);
string findClass(string sub, map<string, set<string> > *MAP);
vector<string> IntersectVec(vector<string> v1, vector<string> v2);
vector<string> IntersectVec2(vector<string> v1, vector<string> v2, int type);
set<string> IntersectSet(set<string> v1, set<string> v2);
set<string> UnionSet(set<string> v1, set<string> v2);
bool contains(string g, set<string> *myset);
bool containsVec(string g, vector<string> *vec);
bool vecContains(int index, vector<int> vec);
bool allNA(vector<double> vec);
bool NAexist(vector<double> vec);
void UniqueVec(vector<string> &v);
set<string> getNonOverlap(vector<string> g1, vector<string> g2);
set<string> findParent(string child, map<string, set<string> >&Map);
int countOcurrence(string str, vector<string> vec);
bool checkNAs(map<string,subjectClass> &subjectMap, int datatype);
int whichpos(string g, vector<string> vec);
vector<vector<double> > removeRow(vector<double> x, vector<vector<double> >&MAT_new);
void checkDuplicate(vector<vector<double> >&MAT);
string getEdge(string str, set<string> *Edges);
void insertPrior(string Priorfile, map<string, subjectClass> &subjectMap);

// statistics functions //
double calMean(vector<double> vec);
double calMedian(vector<double> vec);
double calSdn(vector<double> vec);
double calSd(vector<double> vec);
double calCov(vector<double> x, vector<double> y);
double calPearsonCor(vector<double> x, vector<double> y);
double calEdist(vector<double> v1, vector<double> v2);
double computeL2norm(double dd, vector<double> vec);
double calMeanDik(vector<double> vec);
vector<double> Benjamini_Hochberg(vector<double> pvalue);
void KNNimpute(vector<int> indices,vector<string> subjects, vector<string> genes, vector<vector<double> > MAT,map<set<string>, vector<vector<double> > >&MAP, map<string, subjectClass> &subjectMap, int knnK, int datatype);
vector<double> calweightedEdist( vector<double> vecA, vector<vector<double> > candidateVec);
double imputeMV(int indexNA,vector<double> weightedD, int K, vector<vector<double> > candidateVec);
vector<int> findkminPos(vector<double> v, int k);
vector<int> findkmaxPos(vector<double> v, int k);
vector<double> Ztrans(vector<double> vec, double mean, double sd);

// functions for generating results //
//vector<string> extractInfo(vector<string> GENE,vector<string> SUBJ, set<string> G,map<string, subjectClass>&subjectMap, map<string, set<string> >&MAP, int type);
set<string> identifyGenes(int mink, map<string, set<string> > &subtypeMap, map<string, subjectClass>&subjectMap, int datatype);
void outputData(vector<string> subjects, vector<string> genes, map<string, subjectClass>&subjectMap, int datatype ,string dir);
vector<vector<set<string> > >outputGeneSurv(string dir, bool useZ,bool usePrior,string normby, vector<double> minThres,set<string> *Alledges,vector<string> *subjects, map<string, set<string> > *subtypeMap, map<string, PAM_t> *PAMmap, map<string, subjectClass> *subjectMap, map<string, int> *InteractMap);

//Intermediate functions//
void identifyGenesEDGE(int mink,set<string> genes, map<string, set<string> >&subtypeMap, map<string, subjectClass>&subjectMap, map<string, set<string> > &Map, int code);
vector<string> getSubjects(int datatype, vector<string> genes, map<string, subjectClass> &subjectMap);
vector<double> extractGeneVec(int datatype, string g, vector<string> sub, map<string, subjectClass>&subjectMap);
vector<vector<double> > creatematrix(vector<string> genes, vector<string> subjects, map<string, subjectClass>&subjectMap, int datatype);
vector<string> SetToVec(set<string> SET);
vector<vector<double> > getCandidateG(int pos,vector<double> x_j, vector<vector<double> >MAT);
vector<double> getColvec(int pos, vector<vector<double> > MAT);
void determineCandidate(int max_p,vector<string> genes, vector<vector<double> > MAT, map<set<string>, vector<vector<double> > >&MAP);
vector<double> fillPAMmap(bool useZ, string normby, set<string> *Alledges,vector<string> *sub,map<string, int> *InteractMap, map<string, set<string> > *subtypeMap,map<string, subjectClass> *subjectMap, map<string, PAM_t> &PAMmap, map<string, vector<set<string> > >  *FDneighborMap);
vector<double> fillPAMmap2(bool useZ, string normby, set<string> *Alledges,vector<string> *sub,map<string, int> *InteractMap, map<string, set<string> > *subtypeMap,map<string, subjectClass> *subjectMap, map<string, PAM_t> &PAMmap, map<string, vector<set<string> > > *FDneighborMap, string dir);
double CVKfold(string dir, bool useZ,bool usePrior,string normby,int Kfold,set<string> *Alledges,vector<string> *sub,  int thres_size, map<string, int> *InteractMap,map<string, vector<set<string> > > *NeighborMap, map<string, set<string> > *subtypeMap, map<string, subjectClass>*subjectMap,map<string, PAM_t> *PAMmap);
double calNewXbarik(PAM_t &pam_struct, double thres, int group);
double calNewDik(vector<double> d_ij, double thres, int group);
double getPenalty(double thres, int grp, PAM_t &pam_struct);
double ModifyThres(double thres, double factor, int grp, PAM_t &pam_struct, bool verbose);
set<string> createNeighborMap(string dir,  vector<string> *X_g, vector<string> *Y_g, map<string, set<string> > *BTWmap, map<string, set<string> > *TargetMap,map<string, set<string> > *WITHINmap, map<string, vector<set<string> > >&FDNeighborMap);
vector<vector<string> > GetNeighbors(string str, map<string, vector<set<string> > >&Map);
vector<vector<double> > GetNeighborDik(vector<string> vec, map<string, PAMClass> &pammap);
void StandardizeFeatures(string dir,vector<string> sub, vector<string> genes,int datatype,map<string, subjectClass> &subjectMap); 
double checkCor(string tagA, string tagB, vector<string> sub, map<string, subjectClass>&subjectMap);
double calNumSurvived(double thres, vector<double> maxthres, set<string> *Alledges,map<string, set<string> > *subtypeMap, map<string, PAM_t>*PAMmap);

// functions for Map //
string maximumG(map<string, set<string> > &map);
set<string> retrieveKeys(map<string, subjectClass> &MAP);
void copyMap(map<string, subjectClass> &MapA, map<string, subjectClass> &MapB);
vector<string> getKey(map<string, set<string> > *Map); 
string getTag(string str, map<string, string > *Map); 
void ReverseMap(map<string, set<string> > &Map1,map<string, set<string> >&Map2);
void CleanupMap(vector<string> subjects,map<string, set<string> > &subtypeMap);


// functions for reading in data //
bool readUserInput(string file ,string dir, Input_t &UserInput);
set<string> insertSubtypeMap(string subtypeFile, map<string, set<string> >&subtypeMap);
vector<string> readFile(string myFile,  map<string, set<string> >&subtypeMap,map<string, subjectClass> &subjectMap,Input_t &UserInput, set<string>sub, int dataType);
set<string> readbetweenNetwork(string file,map<string,set<string> >&BTWMap,map<string, int> &InteractMap,vector<string> target_g, vector<string> X_g, set<string> sub, map<string, subjectClass>&subjectMap);
set<string> readwithinNetwork(string file, map<string,set<string> >&WITHINMap,map<string, int> &InteractMap, vector<string> X_g, set<string> sub, map<string, subjectClass>&subjectMap);
set<string> readPATHWAY(string file, map<string,set<string> >&PATHmap, map<string, string>&PATHwayAnnot);

// for Network Enrichment //
void  CreateOutput(string directory,string lab,vector<string > PATHlist, int dir, set<string> *geneLIST, set<string> *bgLIST, map<string, set<string> > *PATHEdgemap,map<string, set<string> > *PATHmap, map<string, string>*PATHwayAnnot,int mink);
vector<double> hypergeo(set<string> *vec, set<string> *DE, set<string> *bglist);
vector<string>createEdgeOrientedPathwayMap(set<string> *features,set<string> *Alledges,map<string, set<string> > *PATHmap,map<string, set<string> > *PATHmapEdge, double prop,int minbg);
set<string> createMEM(set<string> *SET, set<string> *ALLedges);
void printBGlist(string dir, set<string> edges);
#endif /* GLOBALS_HPP */
