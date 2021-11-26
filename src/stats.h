#ifndef STATS_H
#define STATS_H

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <map>
#include "options.h"
#include <iostream>
#include <fstream>
#include <vector>
#include "bed.h"

#define MAX_SUPPORTING_READS 100

class Stats {
public:
    Stats(Options* opt);
    ~Stats();
    void addRead(bam1_t* b);
    void addMolecule(unsigned int supportingReads, bool PE);
    void addCluster(bool hasMultiMolecule);
    void print();
    long getMappedBases();
    long getMappedReads();
    long getReads();
    long getBases();
    void reportJSON(ofstream& ofs);
    double getMappingRate();
    double getDupRate();
    long getMolecules() { return mMolecule; }
    double getMismatchRate();
    void makeGenomeDepthBuf();
    void makeBedStats(Bed* other = NULL);
    void statDepth(int tid, int start, int len);
    void setPostStats(bool flag);
    void addSSCS();
    void addDCS();
    // position unique variation description, similar to "Variant Description String" in VarDictJava.
    void addDupVariety(unsigned int tid, const string& upv);

public:    
	static string list2string(double* list, int size);
    static string list2string(double* list, int size, long* coords);
    static string list2string(long* list, int size);

public:
    long mReadWithMismatches;
    long mCluster;
    long mMultiMoleculeCluster;
    long mMolecule;
    long mMoleculeSE;
    long mMoleculePE;
    long* mSupportingHistgram;
	Options* mOptions;
	long mBase;
	long mBaseMismatches;
	long mBaseUnmapped;
	long mRead;
	long mReadUnmapped;
    long uncountedSupportingReads;
    bool mIsPostStats;
    long mSSCSNum;
    long mDCSNum;
    vector<vector<long>> mGenomeDepth;
    // mutation occurs in multiple different duplicated group statistic, two levels' map make the output sorted
    map<unsigned int, map<string, int>> varDupVariety;
    Bed* mBedStats;
};

#endif