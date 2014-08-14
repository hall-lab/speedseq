#ifndef __INTERVAL_HH__
#define __INTERVAL_HH__

// C/C++ includes
#include <iostream> 
#include <fstream> 
#include <sstream> 
using namespace std; 

class Interval
{
private:
  string _name1,_name2;
  short _index;
  int _s1,_e1;
  int _s2,_e2;
  short _u1,_u2;
  Interval *_next;

public:
  Interval(string file);

private:
  Interval();

public:
  inline string name1()      { return _name1; }
  inline int    start1()     { return _s1; }
  inline int    end1()       { return _e1; }
  inline int    unique1()    { return _u1; }
  inline string name2()      { return _name2; }
  inline int    start2()     { return _s2; }
  inline int    end2()       { return _e2; }
  inline int    unique2()    { return _u2; }
  inline int    otherIndex() { return _index; }
  inline Interval *next()    { return _next; }

public:
  void setUnique1(int val);
  void setUnique2(int val);
  void setOtherIndex(short val);
  int getEquivalent(int pos);
  int count();

private:
  inline void setNext(Interval *inter) { _next = inter; }
  bool parseRecord(ifstream &in);
};

#endif
