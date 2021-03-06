// written and debugged by Joel Graber, Senior Research Associate at the
// Center for Advanced Biotechnology, Boston University
// All rights reserved by The Trustees of Boston University, 1999.
//
#ifndef __LOGFILE_H
#define __LOGFILE_H
//
#ifndef _STRING_
#include <string>
#endif
#ifndef _FUNCTIONAL_
#include <functional>
#endif
#include <vector>
using std::string; using std::unary_function; using std::binary_function; using std::vector;
typedef vector<string> strVec;
//////////////////////////////////////////////////////////////////////////////////////////////
class LogSetup:public binary_function<int,char**,bool> {
 public:
  bool operator()(int argc,char *argv[]) const;
};
//////////////////////////////////////////////////////////////////////////////////////////////
class InitLog:public binary_function<string, bool, int> {
 public:
  int operator()(const string& fn,bool overwrite=true) const;
};
//////////////////////////////////////////////////////////////////////////////////////////////
class InitErrLog:public binary_function<string, bool, int> {
 public:
  int operator()(const string& fn,bool overwrite=true) const;
};
//////////////////////////////////////////////////////////////////////////////////////////////
class WriteToLogTS:public unary_function<string,int> {
 public:
  int operator()(const string& s) const;
}; //
//////////////////////////////////////////////////////////////////////////////////////////////
class CloseLogTS:public unary_function<string,int> {
 public:
  int operator()(const string& s) const;
}; //
//////////////////////////////////////////////////////////////////////////////////////////////
class WriteToLog:public unary_function<string,int> {
 public:
  int operator()(const string& s) const;
}; //
//////////////////////////////////////////////////////////////////////////////////////////////
class WriteIntToLog:public binary_function<string,int,int> {
 public:
  int operator()(const string& s,int i) const;
}; //
//////////////////////////////////////////////////////////////////////////////////////////////
class WriteFloatToLog:public binary_function<string,float,int> {
 public:
  int operator()(const string& s,float f) const;
}; //
//////////////////////////////////////////////////////////////////////////////////////////////
class WriteToErrTS:public unary_function<string,int> {
 public:
  int operator()(const string& s) const;
}; //
//////////////////////////////////////////////////////////////////////////////////////////////
class CloseErrTS:public unary_function<string,int> {
 public:
  int operator()(const string& s) const;
}; //
//////////////////////////////////////////////////////////////////////////////////////////////
class WriteToErr:public unary_function<string,int> {
 public:
  int operator()(const string& s) const;
}; //
//////////////////////////////////////////////////////////////////////////////////////////////
class WriteIntToErr:public binary_function<string,int,int> {
 public:
  int operator()(const string& s,int i) const;
}; //
//////////////////////////////////////////////////////////////////////////////////////////////
class WriteFloatToErr:public binary_function<string,float,int> {
 public:
  int operator()(const string& s,float f) const;
}; //
//////////////////////////////////////////////////////////////////////////////////////////////
class WriteToLogs:public unary_function<string,int> {
 public:
  int operator()(const string& s) const;
}; //
//////////////////////////////////////////////////////////////////////////////////////////////
class WriteIntToLogs:public binary_function<string,int,int> {
 public:
  int operator()(const string& s,int i) const;
}; //
//////////////////////////////////////////////////////////////////////////////////////////////
class WriteFloatToLogs:public binary_function<string,float,int> {
 public:
  int operator()(const string& s,float f) const;
}; //
//////////////////////////////////////////////////////////////////////////////////////////////
class WriteRunLog:public unary_function<strVec,bool> {
 public:
  bool operator()(const strVec& a) const;
};
//////////////////////////////////////////////////////////////////////////////////////////////
bool WroteError();
void LogsOff();
//////////////////////////////////////////////////////////////////////////////////////////////
#endif  // __LOGFILE_H
