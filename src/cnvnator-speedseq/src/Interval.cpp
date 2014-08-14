#include "Interval.hh"

Interval::Interval() : _name1(""),_s1(0),_e1(0),_u1(0),
		       _name2(""),_s2(0),_e2(0),_u2(0),
		       _index(-1),_next(NULL)
{}

Interval::Interval(string file) : _name1(""),_s1(0),_e1(0),_u1(0),
				  _name2(""),_s2(0),_e2(0),_u2(0),
				  _index(-1),_next(NULL)
{
  ifstream in(file.c_str());
  if (!in.good()) {
    cerr<<"Can't open file '"<<file<<"'."<<endl;
    return;
  }
  if (!parseRecord(in)) {
    in.close();
    return;
  }
  Interval *curr = this;
  while (!in.eof()) {
    Interval *next = new Interval();
    if (next->parseRecord(in)) {
      curr->setNext(next);
      curr = next;
    } else {
      delete next;
      break;
    }
  }
  in.close();
  return;  
}

bool Interval::parseRecord(ifstream &in)
{
  if (in.eof()) return false;
  in>>_name1;
  if (in.eof()) return false;
  in>>_s1;
  if (in.eof()) return false;
  in>>_e1;
  if (in.eof()) return false;
  in>>_name2;
  if (in.eof()) return false;
  in>>_s2;
  if (in.eof()) return false;
  in>>_e2;
  return true;
}

int Interval::getEquivalent(int pos)
{
  if (pos < _s1 || pos > _e1) return -1;
  if (_e2 > _s2) return (pos - _s1) + _s2;
  else           return _s2 - (pos - _s1);
}

int Interval::count()
{
  Interval *i = this;
  int ret = 0;
  while (i) { i = i->next(); ret++; }
  return ret;
}

void Interval::setUnique1(int val)
{
  if (val < 0 || val > 30000) val = 30000;
  _u1 = val;
}

void Interval::setUnique2(int val)
{
  if (val < 0 || val > 30000) val = 30000;
  _u2 = val;
}

void Interval::setOtherIndex(short val) { _index = val; }
