/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.cpp
 * Author: Hiromi WL Koh
 *
 * Last modified on 13 June, 2018, 4:00 PM
 */

#include "globals.hpp"
#include "geneClass.hpp"
#include "subjectClass.hpp"

#include <algorithm>
#include <cstdlib>
#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <set>
#include <ctime>
#include <cmath>

using namespace std;
/*
 * 
 */
int main(int argc, char** argv) {
    
    if(argc <2){
        cerr <<"\nUSAGE: ./iOmicsPASS <<inputFile> <<directory>>\n";
        return 0;
    }

    clock_t begin = clock();
    clock_t time;
    double elapse, elapse2;

    string inputF = argv[1];	
    string inputD = argv[2];

    Input_t *UserInput = new Input_t();    
    map<string, subjectClass> subjectMap;
    map<string, PAM_t> PAMmap;
    map<string, set<string> > WithinMap;
    map<string, set<string> > BTWMap, PATHmap;
    map<string, set<string> > subtypeMap;
    map<string, vector<set<string> > > FirstDNeighborMap;
    map<string, int> InteractMap;
    map<string, string> PATHwayAnnot;

    vector<string> Z_g,Y_g, X_g;
    set<string> BTW_edge, WITHIN_edge, pathways;
    vector<string> GENE, gene_interest;
    
    cerr << "\n\nInitiating iOmicsPASS...\n" <<endl; 
    readUserInput(inputF, inputD, *UserInput);  
    
    cerr<<"\nApplying the following filter (specified by user) across subtypes:"<<endl;
    cerr<<"Minimum number of non-missing observations for each feature within each subtype: "<< UserInput->minobs<<endl;
    cerr<<"Minimum proportion of non-missing observations for each feature within each subtype: "<<UserInput->minprop<<endl; 

    cerr << "\nReading in the subtype information file....." <<endl;
    set<string> subjects = insertSubtypeMap(UserInput->subtypefile, subtypeMap);    
    vector<string> grplab = getKey(&subtypeMap);  

    cerr<<"\nReading in the data file......"<<endl;
    X_g = readFile(UserInput->dataX,  subtypeMap,subjectMap,*UserInput,subjects, 3);
    cerr<<"\nThere are "<< X_g.size()<<" features in the data X."<<endl;

    if(UserInput->analyzeY){
      Y_g = readFile(UserInput->dataY, subtypeMap, subjectMap,*UserInput,subjects, 2);
      cerr<<"\nThere are "<<Y_g.size()<<" features in data Y."<<endl;
    }

    if(UserInput->analyzeZ){
       Z_g = readFile(UserInput->dataZ, subtypeMap,subjectMap,*UserInput,subjects, 1);
       cerr<<"\nThere are "<<Z_g.size()<<" features in data Z."<<endl;
    }

    cerr << "\nReading in the network file(s)......"<<endl;
    //if Z file is provided, we require the genes to be present in both Y and Z.
    vector<string> target_g;
    if(!UserInput->analyzeZ & UserInput->analyzeY) target_g = Y_g;
    if(UserInput->analyzeZ){
	if(UserInput->analyzeY){
		if(UserInput->normalize_Zby=="Y")  target_g = IntersectVec2(Y_g, Z_g,2);
		if(UserInput->normalize_Zby=="X")  X_g = IntersectVec2(X_g, Z_g,3);
	}
	if(!UserInput->analyzeY) X_g = IntersectVec2(X_g, Z_g,3);
    }
    //if only data X is supplemented, then only within network is needed. If data Y is also supplemented, then between network is needed.
    WITHIN_edge = readwithinNetwork(UserInput->withinNet, WithinMap,InteractMap, X_g, subjects, subjectMap); 
    if(UserInput->analyzeY) BTW_edge = readbetweenNetwork(UserInput->betweenNet, BTWMap,InteractMap, target_g, X_g, subjects, subjectMap);
    set<string> Alledges = UnionSet(WITHIN_edge,BTW_edge);
    printBGlist(inputD,Alledges);

    if(UserInput->enrichment){
    	cerr<<"\nReading in the Pathway module file....."<<endl;
    	pathways = readPATHWAY(UserInput->pathwayfile, PATHmap, PATHwayAnnot);
    }
    if(UserInput->prior_input){
      cerr << "\nPrior is supplemented by user and will be used in classifying samples....."<<endl;
      insertPrior(UserInput->priorfile, subjectMap);
    }
    set<string> allS = retrieveKeys(subjectMap);
    cerr<<"There are "<<allS.size()<<" subjects in the subjectMap"<<endl;

    cerr <<"\nAll data files have been successfully read into iOmicsPASS!\n\nProceeding to next step...\n";
    vector<string> sub_Z, sub_Y, sub_X;
    if(UserInput->analyzeZ) sub_Z = getSubjects(1, Z_g, subjectMap);
    if(UserInput->analyzeY) sub_Y = getSubjects(2, Y_g, subjectMap);
    sub_X = getSubjects(3, X_g, subjectMap);

    cerr<<"There are "<<sub_X.size()<<" subjects with subtype information and "<< X_g.size()<<" features after filtering in data X"<<endl;
    if(UserInput->analyzeY)cerr<<"There are "<<sub_Y.size()<<" subjects with subtype information and "<<Y_g.size()<< " features after filtering in data Y"<<endl;
    if(UserInput->analyzeZ) cerr<<"There are "<< sub_Z.size()<< " subjects with subtype information and "<< Z_g.size()<<" features after filtering in data Z"<<endl;

    bool impute=UserInput->knnimpute;
    bool checkD, checkR, checkP;
    if(impute){
       if(!UserInput->analyzeZ) checkD = false;
       if(UserInput->analyzeZ) checkD=checkNAs(subjectMap,1);
       if(!UserInput->analyzeY) checkR= false;
       if(UserInput->analyzeY) checkR=checkNAs(subjectMap,2);
       checkP=checkNAs(subjectMap,3);
       string new_lab;
       vector<string> lab;
       if(checkD) lab.push_back("Z");
       if(checkR) lab.push_back("Y");
       if(checkP) lab.push_back("X");
       if(lab.size()==0) cerr<<"\nNo missing values detected, skipping KNN-imputation step.\n"<<endl;
       else{
         for(int j=0; j<lab.size(); j++){
  	    if(j==0) new_lab = lab.at(j);
 	    else new_lab+= " and " + lab.at(j);
         }
         cerr<<"\nKNN imputation will be carried out on data "<<new_lab<<".\n"<<endl;
         int knnK = UserInput->knnK;
	 int blocksize = UserInput->blocksize;
	 map<set<string>, vector<vector<double> > > CandidateSET_d, CandidateSET_r, CandidateSET_p;

         if(checkD){
		cerr<<"Imputing data Z...\n";
		vector<vector<double> > MAT_D = creatematrix(Z_g, sub_Z, subjectMap,1);
		determineCandidate(blocksize,Z_g, MAT_D, CandidateSET_d);
		time = clock();
    		elapse2 = double(time - begin)/(double)CLOCKS_PER_SEC;
		KNNimpute(UserInput->cid_Z,sub_Z, Z_g, MAT_D,CandidateSET_d, subjectMap,knnK,1 );
		outputData(sub_Z, Z_g, subjectMap,1, inputD);
		time = clock();
    		elapse = double(time - begin)/(double)CLOCKS_PER_SEC;
    		cerr<<"\nTook : "<<(elapse-elapse2)<<" seconds.\n"<<endl; 

	 }
         if(checkR){
		cerr<<"Imputing data Y...\n";
	  	vector<vector<double> > MAT_R = creatematrix( Y_g,sub_Y, subjectMap, 2);
		determineCandidate(blocksize,Y_g, MAT_R, CandidateSET_r);
		time = clock();
    		elapse2 = double(time - begin)/(double)CLOCKS_PER_SEC;
		KNNimpute(UserInput->cid_Y,sub_Y, Y_g, MAT_R, CandidateSET_r,subjectMap, knnK,2);
		outputData(sub_Y, Y_g, subjectMap,2,inputD);
		time = clock();
     		elapse = double(time - begin)/(double)CLOCKS_PER_SEC;
    		cerr<<"\nTook : "<<(elapse-elapse2)<<" seconds.\n"<<endl; 

	 }
         if(checkP){
		cerr<<"Imputing data X...\n";
		vector<vector<double> > MAT_P = creatematrix( X_g,sub_X,subjectMap,3);
		determineCandidate(blocksize, X_g, MAT_P, CandidateSET_p);
		time = clock();
    		elapse2 = double(time - begin)/(double)CLOCKS_PER_SEC;
		KNNimpute(UserInput->cid_X,sub_X,X_g, MAT_P, CandidateSET_p, subjectMap, knnK,3);
		outputData(sub_X, X_g, subjectMap,3, inputD);	
		time = clock();
    		elapse = double(time - begin)/(double)CLOCKS_PER_SEC;
    		cerr<<"\nTook : "<<(elapse-elapse2)<<" seconds.\n"<<endl; 

         }
       }
    }   
    vector<string> commonSub;
    if(!UserInput->analyzeY) commonSub =sub_X;
    if(UserInput->analyzeY) commonSub = IntersectVec(sub_Y, sub_X);
    if(UserInput->analyzeZ) commonSub = IntersectVec(commonSub,sub_Z);

    // standardize measurements in the map
    if(UserInput->analyzeZ & UserInput->ztrans_Z) StandardizeFeatures(inputD, commonSub,Z_g,1,subjectMap);
    if(UserInput->ztrans_Y & UserInput->analyzeY) StandardizeFeatures(inputD, commonSub, Y_g,2, subjectMap);
    if(UserInput->ztrans_X) StandardizeFeatures(inputD, commonSub, X_g,3, subjectMap);

    cerr<<"\n\nNumber of common subjects across all data is "<<commonSub.size()<<endl;
 
    CleanupMap(commonSub, subtypeMap);
    map<string, set<string> > Targets2BTWmap;
    if(UserInput->analyzeY) ReverseMap(BTWMap, Targets2BTWmap);
    set<string> features = createNeighborMap(inputD, &X_g, &target_g, &BTWMap,&Targets2BTWmap, &WithinMap, FirstDNeighborMap);
    cerr<<"Carrying out shrunken centroid module on the co-modulating edges..."<<endl;

    bool useZ = UserInput->analyzeZ;
    vector<double> thresMax = fillPAMmap2(useZ, UserInput->normalize_Zby, &Alledges,&commonSub, &InteractMap, &subtypeMap,&subjectMap, PAMmap, &FirstDNeighborMap, inputD);
    cerr<<"\nMaximum Dij (average) in the data is : "<< calMean(thresMax) <<"\n"<<endl;
    cerr<<"Class-specific thresholds are:";
    for(int i=0; i<thresMax.size(); i++) cerr<<"\t"<<thresMax.at(i);
    cerr<<endl;

    double minTHRES, minThres;
    minTHRES = UserInput->minthres;
    if(UserInput->crossV){
	cerr<<"\nCarrying out cross-validation...(This may take a while)\n"<<endl;
    	minThres = CVKfold(inputD, useZ,UserInput->prior_input,UserInput->normalize_Zby,UserInput->numFold, &Alledges,&commonSub,  30, &InteractMap, &FirstDNeighborMap, &subtypeMap, &subjectMap, &PAMmap);
    }
    if(!UserInput->crossV) cerr<<"\nSkipping cross-validation and using user-specified threshold: "<<minTHRES<<"\n"<<endl;
    if(minTHRES==NA_VAL) minTHRES = minThres;
    vector<double> MODthres = thresMax;
    for(int i=0; i<thresMax.size(); i++){
        MODthres.at(i) = (minTHRES/calMean(thresMax))*thresMax.at(i);
    }
    vector<vector<set<string> > > geneSurv = outputGeneSurv(inputD, useZ,UserInput->prior_input,UserInput->normalize_Zby, MODthres,&Alledges, &commonSub, &subtypeMap,  &PAMmap, &subjectMap, &InteractMap);

    //Network-centric Enrichment
    if(!UserInput->enrichment) cerr<<"\nSkipping sub-network enrichment module.\n"<<endl;
    if(UserInput->enrichment){
	   cerr<<"\nCarrying out sub-network enrichment on the predictive edges for each subtype..."<<endl;
    	   map<string, set<string> > PATHmapEdge;
           vector<string> PATHs = createEdgeOrientedPathwayMap(&features,&Alledges, &PATHmap, &PATHmapEdge,UserInput->bgProp, UserInput->minbgSize);
           set<string> EdgeSurv_grp;

           for(int j=0; j<2;j++){
      		for(int i=0; i<grplab.size(); i++){
         	    EdgeSurv_grp = (geneSurv.at(j)).at(i);
         	    CreateOutput(inputD,grplab.at(i), PATHs, j, &EdgeSurv_grp, &Alledges, &PATHmapEdge, &PATHmap,&PATHwayAnnot, UserInput->minsigSize);
      	    	}
    	   }
    }
 

    cerr <<"\n\nThe program has finished running!\n" ;
    clock_t end = clock();
    double elapse_sec = double(end - begin)/(double)CLOCKS_PER_SEC;
    string timetype;
    double new_t, tmpt;
    tmpt = elapse_sec/60.0;
    double tt = round( tmpt*10)/10;
    if(tt>240){
	timetype = " hours";
	new_t = round((tmpt/60.0)*10)/10;
    }
    if(tt<=240){
	timetype =" minutes";
	new_t = tt;
    }
    cerr<<"\nTime to completion of iOmicsPASS is : "<<new_t<<timetype<<endl;
    cerr<<"----------------------------------------------------------------"<<endl;
    
    return 0;
}

