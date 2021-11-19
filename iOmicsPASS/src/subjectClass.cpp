/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   subjectClass.cpp
 * Author: Hiromi WL Koh
 * 
 * Last modified on 18 February, 2020, 11:00 AM
 *
 */


#include "globals.hpp"
#include "geneClass.hpp"

#include <algorithm>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

subjectClass::subjectClass() {
}

subjectClass::subjectClass(const subjectClass& orig) {
}

subjectClass::~subjectClass() {
}

subjectClass::subjectClass(string id){
    subject = id;
}

set<string> subjectClass::getCommonG(int datatype1, int datatype2){
     set<string> setA,setB;
     if(datatype1==1) setA = Z;
     if(datatype1==2) setA = Y;
     if(datatype1==3) setA = X;
     if(datatype2==1) setB = Z;
     if(datatype2==2) setB = Y;
     if(datatype2==3) setB = X;
     set<string> Iset = IntersectSet(setA,setB);
     return Iset;
}

void subjectClass::removeNAgenes(string g, int datatype){
   set<string>::iterator it;
   if(datatype==1){
	for(it =Z_missing.begin(); it!=Z_missing.end(); it++) if(*it==g) Z_missing.erase(it);
   }
   if(datatype==2){
	for(it = Y_missing.begin(); it!=Y_missing.end(); it++) if(*it==g) Y_missing.erase(it);
   }
   if(datatype==3){
	for(it = X_missing.begin(); it!=X_missing.end(); it++) if(*it==g) X_missing.erase(it);
   }
}

void subjectClass::insertNAgenes(string g, int datatype){
  
     if(datatype==1) Z_missing.insert(g);
     if(datatype==2) Y_missing.insert(g);
     if(datatype==3) X_missing.insert(g);
}

set<string> subjectClass::getMissingG(int datatype){
     set<string> ret;
     if(datatype==1) ret = Z_missing;
     if(datatype==2) ret = Y_missing;
     if(datatype==3) ret = X_missing;
     return ret;
}


void subjectClass::insertMap(string g, double val ,int dataType){
    
    geneClass *curG = NULL;
    map<string,geneClass>::iterator m_iter, m_iter2;
    m_iter = geneMap.find(g);
    
    if(m_iter==geneMap.end()){
        curG = new geneClass(g);
        geneMap[g]= *curG;
        delete curG;
        m_iter2 = geneMap.find(g);
        if(dataType==1 & val!= NA_VAL) {
            m_iter2->second.addZ(val);
            Z.insert(g);
	    Z_exist = true;
        }
        if(dataType==2 & val!= NA_VAL) {
            m_iter2->second.addY(val);
            Y.insert(g);
	    Y_exist = true;
        }
        if(dataType==3 & val!= NA_VAL) {
            m_iter2->second.addX(val);
            X.insert(g);
	    X_exist = true;
        }
    }
    else{
        if(dataType==1 & val!= NA_VAL) {
            m_iter->second.addZ(val);
            Z.insert(g);
	    Z_exist = true;
        }
        if(dataType==2 & val!= NA_VAL) {
            m_iter->second.addY(val);
            Y.insert(g);
	    Y_exist = true;
        }
        if(dataType==3 & val!= NA_VAL) {
            m_iter->second.addX(val);
            X.insert(g);
	    X_exist = true;
        }
    }
}
void subjectClass::insertPrior(string grp, double prob){
          Group.push_back(grp);
	  probs.push_back(prob);
	  if(prob!= NA_VAL) prior_exist = true;
}
double subjectClass::getPrior(string grp){
	double ret=1.0;
	bool found = false;
        int index_g;
	for(int i=0; i<Group.size();i++){
		string currG = Group.at(i);
		if(currG==grp){
		   index_g = i;	
		   found = true;
		}
		if(found) continue;
	}
	if(found) ret = probs.at(index_g);
	return ret;
}
set<string> subjectClass::getSet(int dataType){
    set<string> ret;
    if(dataType==1) ret = Z;
    if(dataType==2) ret = Y;
    if(dataType==3) ret = X;
    return ret;
}

bool subjectClass::existSet(int dataType){
    bool ret;
    if(dataType==1) ret = Z_exist;
    if(dataType==2) ret = Y_exist;
    if(dataType==3) ret = X_exist;
    return ret;
}



