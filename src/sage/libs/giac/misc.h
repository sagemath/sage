#ifndef GIACPYMISC_H
#define GIACPYMISC_H
#include <giac/giac.h>

#include <fstream>

using namespace std;
using namespace giac;

inline void  ressetctrl_c(void){
giac::ctrl_c=false;
}

inline int  testctrl_c(void){
  if(giac::ctrl_c){
    return 1;}
  else{
    return 0;
  }
}

/*
int  giacgencmp(gen & a, gen & b, const context * context_ptr){
  if(a==b){return 0;}
  else{ if (is_strictly_greater(a,b,context_ptr)){return 1;}
    else{return -1;}
  }
}
*/

int  giacgenrichcmp(gen & a, gen & b, int op, const context * context_ptr){
  int rep=0;//be careful with undef results
  switch ( op )
    {
    case 0 : //<
       if (is_strictly_greater(b,a,context_ptr)) {rep = 1;}
       break;
    case 1 : //<=
       if (is_greater(b,a,context_ptr)) {rep = 1;}
       break;
    case 2 ://==
       if (operator_equal(b,a,context_ptr)) {rep = 1;}
       break;
    case 3 ://!=
       if (!operator_equal(b,a,context_ptr)) {rep = 1;}
       break;
    case 4 ://>
       if (is_strictly_greater(a,b,context_ptr)) {rep = 1;}
       break;//>=
     case 5 :
       if (is_greater(a,b,context_ptr)) {rep = 1;}
       break;
     default :
       rep=0;
    };
  return rep;
}

gen  giacmul(gen & a, gen & b, const context * context_ptr){
  return eval(a*b,context_ptr);
}

gen  giacdiv(gen & a, gen & b, const context * context_ptr){
  return eval(a/b,context_ptr);
}

gen  giacmod(gen & a, gen & b, const context * context_ptr){
  if(b != 0)
     return eval(a * makemod(1,b),context_ptr);
  else
    return eval(makemod(a,b),context_ptr);
}

int htmlbrowserhelp(char * s){
  if (system_browser_command(s))
    { return 0;}
  else
    {return 1;}
}

string browser_help(const giac::gen & g, int language){
    giac::gen f(g);
    string s;
    giac::html_help_init("aide_cas",language,true,true);
    if (f.type==giac::_SYMB)
      f=f._SYMBptr->sommet;
    if (f.type==giac::_FUNC)
      s=f._FUNCptr->ptr()->s;
    giac::html_vtt=giac::html_help(giac::html_mtt,s);
     if (!giac::html_vtt.empty()){
       //giac::system_browser_command(giac::html_vtt.front());
       return giac::html_vtt.front();
     }
     else{
       return "";
     }
}


void archivegen( const string filename, const gen & g, const context * context_ptr){
  ofstream of(filename.c_str());
  giac::archive(of,g,context_ptr);
  of.close();
}

gen unarchivegen( const string filename, const context * context_ptr){
  ifstream f(filename.c_str());
  gen g=giac::unarchive(f,context_ptr);
  f.close();
  return g;
}


#endif

