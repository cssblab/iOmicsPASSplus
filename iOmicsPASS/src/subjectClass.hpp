/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   subjectClass.hpp
 * Author: Hiromi WL Koh
 *
 * Last modified on 18 February, 2020, 11:00 AM
 */

#ifndef SUBJECTCLASS_HPP
#define SUBJECTCLASS_HPP

#include "globals.hpp"
#include "geneClass.hpp"

#include <vector>
#include <map>
#include <set>
#include <string>
#include <iostream>
#include <fstream>

using namespace std;


class subjectClass {
    
private:
    string subject;
    bool Z_exist=false;
    bool Y_exist=false;
    bool X_exist=false;
    bool prior_exist = false;
    set<string> Z;
    set<string> Y;
    set<string> X;
    set<string> Z_missing, Y_missing, X_missing;
    vector<string> Group;
    vector<double> probs;

public:
    subjectClass();
    subjectClass(const subjectClass& orig);
    virtual ~subjectClass();
    
    subjectClass(string id);
    void removeNAgenes(string g, int datatype);
    void insertMap(string g, double val ,int dataType);
    void insertNAgenes(string g, int datatype);
    void insertPrior(string grp, double prob);
    double getPrior(string grp);
    set<string> getMissingG(int datatype);
    set<string> getSet(int dataType);
    set<string> getCommonG(int datatype1, int datatype2);
    bool existSet(int type);
    map<string, geneClass> geneMap;
    bool checkPrior(){return prior_exist;}
    vector<string> getGrp(){return Group;}
    vector<double> getPriorVec(){return probs;}
};

#endif /* SUBJECTCLASS_HPP */

