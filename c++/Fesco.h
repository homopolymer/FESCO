#ifndef FESCO_H
#define FESCO_H

#include <ctime>
#include <string>
#include <iomanip>
#include <iostream>
#include <vector>
#include <unordered_map>
using namespace std;

void runtime_message(const string& msg, bool verbose);


int fesco_inverted_index_build(int argc, char* argv[]);
int fesco_inverted_index_stats(int argc, char* argv[]);
int fesco_inverted_index_sketch(int argc, char* argv[]);
int fesco_inverted_index_cossim(int argc, char* argv[]);

void fesco_save_sketches(const string& filename,const unordered_map<string,vector<string>>& read_sketches);
void fesco_load_sketches(const string& filename,unordered_map<string,vector<string>>& read_sketches);
#endif
