#include "Fesco.h"
#include <string>
#include <iostream>
using namespace std;

void runtime_message(const string& msg, bool verbose)
{
    time_t result = time(nullptr);
    char mbstr[128];
    strftime(mbstr,sizeof(mbstr),"%c %Z", localtime(&result));

    string msg2 = "[" + string(mbstr) + "]: " + msg + "\n";

    if (verbose) cerr << msg2;
}

void help()
{
    cerr << "Fast and error-tolerant long sequence comparison" << endl;
    cerr << "" << endl;
    cerr << "Usage: fesco <command> <options>" << endl;
    cerr << "" << endl;
    cerr << "Commands: build    Build inverted index for sequences" << endl;
    cerr << "          sketch   Compute sketchs for sequences" << endl;
    cerr << "          cossim   Compute cosine similarity for sequences" << endl;
    cerr << "          stats    Statistical summary of an inverted index" << endl;
    cerr << "" << endl;
    cerr << "--help,-h          Print this message" << endl;
}


int main(int argc, char* argv[])
{
    if (argc==1 || string(argv[1])=="-h" || string(argv[1])=="--help"){
        help();
        return 0;
    }else if (string(argv[1])=="build"){
        return fesco_inverted_index_build(--argc,++argv);
    }else if (string(argv[1])=="stats"){
        return fesco_inverted_index_stats(--argc,++argv);
    }else if (string(argv[1])=="sketch"){
        return fesco_inverted_index_sketch(--argc,++argv);
    }else if (string(argv[1])=="cossim"){
        return fesco_inverted_index_cossim(--argc,++argv);
    }else{
        cerr << "Unknown command: " << argv[1] << endl;
        return -1;
    }
}
