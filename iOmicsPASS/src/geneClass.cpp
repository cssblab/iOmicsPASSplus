/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   geneClass.cpp
 * Author: Hiromi WL Koh
 * 
 * Last modified on 13 June, 2018, 4:00 PM
 */

#include "globals.hpp"

#include <algorithm>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

geneClass::geneClass() {
}

geneClass::geneClass(const geneClass& orig) {
}

geneClass::~geneClass() {
}

geneClass::geneClass(string id){
    gene = id;
}

void geneClass::addX(double p){
    x = p;
    X_bool = true;
}

void geneClass::addY(double r){
    y = r;
    Y_bool = true;
}

void geneClass::addZ(double d){
    z = d;
    Z_bool = true;
}

bool geneClass::ratioAvail(int datatype){
  bool ret = false;
	if(datatype==1) ret = Z_bool;
	if(datatype==2) ret = Y_bool;
	if(datatype==3) ret = X_bool;
	return ret;
}

double geneClass::getReadings(int datatype){
    double ret;
    if(datatype==1) {
	if(Z_bool) ret = z;
	else ret=NA_VAL;
    } 
    if(datatype==2) {
	if(Y_bool) ret = y;
	else ret = NA_VAL;
    }
    if(datatype==3) {
	if(X_bool) ret = x;
	else ret = NA_VAL;
    }
    return ret;
}

