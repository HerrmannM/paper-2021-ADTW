/** Copyright ADTW authors - This is a demonstration application only. **/
#include <iostream>
#include <vector>
#include <string>
#include <filesystem>
using namespace std;
namespace fs = std::filesystem;


// --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
// Argument checking/parsing
// --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
struct Args{
  fs::path ucr_path;
  string dsname;
  double penalty;
};

[[noreturn]] void arg_error(){
  cerr << "Argument error" << endl;
  cerr << "./path/exec <ucr folder> <dataset name> <penalty>" << endl;
  exit(1);
}

Args read_args(int argc, char** argv){
  vector<string> argList(argv, argv+argc);
  Args res;
  int i=1;
  // ---
  if(i<argc){ res.ucr_path = fs::path(argList[i]); ++i; } else { arg_error(); }
  // ---
  if(i<argc){ res.dsname = argList[i]; ++i; } else { arg_error(); }
  // ---
  if(i<argc){ res.penalty = stod(argList[i]); ++i; } else { arg_error(); }
  // ---
  return res;
}


// --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
// Main app
// --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

int main(int argc, char** argv) {
  Args args = read_args(argc, argv);
  // --- Load data - bailout if contains missing data, accept unequal length
  // ---

  return 0;
}





