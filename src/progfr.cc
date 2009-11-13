// -*- mode:C++ ; compile-command: "g++ -I.. -fPIC -DPIC -g -c progfr.cc -o progfr.lo && ln -sf progfr.lo progfr.o && gcc -shared progfr.lo -lc  -Wl,-soname -Wl,libprogfr.so.0 -o libprogfr.so.0.0.0 && ln -sf libprogfr.so.0.0.0 libprogfr.so.0 && ln -sf libprogfr.so.0.0.0 libprogfr.so" -*-
#include "first.h"
/*
 *  Copyright (C) 2002 B. Parisse, Inst. Fourier, 38402 St Martin d'Heres
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */
/* Implementation pour xcas d'un langage d'algorithmique en francais
   Les differences par rapport au "gros rouge"de Renee De Graeve
   http://www-fourier.ujf-grenoble.fr/~degraeve/grouge.pdf
   sont les suivantes:
   * il y a un separateur d'instructions ; (point-virgule)
   * les listes sont delimitees par [] et non {}
   * l'affectation se fait par := au lieu de -> (en sens inverse)
   * pour la valeur de fin d'une boucle le mot-clef est jusque
   * pour une boucle indefinie utiliser tantque (pas d'espace)
   * pas d'accent dans resultat pour la valeur de retour d'une fonction
   * la definition d'une fonction se fait par f(x,y,...):={ local a,b,c; ...}
   */
using namespace std;
#include <giac/input_parser.h>
#include <stdexcept>
#include <cmath>
#include <cstdlib>
#include <giac/giac.h>

#ifndef NO_NAMESPACE_GIAC
namespace giac {
#endif // ndef NO_NAMESPACE_GIAC

  string printassialorssinon(const gen & feuille,const string & sommetstr){
    if ( (feuille.type!=_VECT) || (feuille._VECTptr->size()!=3) )
      return "sialorssinon("+feuille.print()+')';
    const_iterateur it=feuille._VECTptr->begin(),itend=feuille._VECTptr->end();
    string res("si ");
    res += sametoequal(*it).print();
    ++it;
    res += " alors ";
    indent_spaces +=2;
    if ((it->type==_SYMB) && (it->_SYMBptr->sommet==at_bloc))
      res += printasinnerbloc(it->_SYMBptr->feuille);
    else
      res += it->print() ;
    indent_spaces -=2;
    res+= " sinon ";
    ++it;
    indent_spaces +=2;
    if ( (it->type==_SYMB) && (it->_SYMBptr->sommet==at_bloc))
      res += printasinnerbloc(it->_SYMBptr->feuille);
    else
      res += it->print() ;
    indent_spaces -=2;
    res += indent()+ "fsi;";
    return res;
  }

  gen _sialorssinon(const gen & g){
    return _ifte(g);
  }
  const string _sialorssinon_s("sialorssinon");
  unary_function_unary __sialorssinon(&_sialorssinon,_sialorssinon_s,&printassialorssinon);
  unary_function_ptr at_sialorssinon (&__sialorssinon,_QUOTE_ARGUMENTS,T_IFTE);

  const string _si_s("si");
  unary_function_unary __si(&_sialorssinon,_si_s,&printassialorssinon);
  unary_function_ptr at_si (&__si,_QUOTE_ARGUMENTS,T_IF);
  
  const string _alors_s("alors");
  unary_function_unary __alors(&_ifte,_alors_s);
  unary_function_ptr at_alors (&__alors,_QUOTE_ARGUMENTS,T_THEN);
  
  const string _sinon_s("sinon");
  unary_function_unary __sinon(&_ifte,_sinon_s);
  unary_function_ptr at_sinon (&__sinon,_QUOTE_ARGUMENTS,T_ELSE);
  
  const string _fsi_s("fsi");
  unary_function_unary __fsi(&_ifte,_fsi_s);
  unary_function_ptr at_fsi (&__fsi,_QUOTE_ARGUMENTS,T_BLOC_END);

  string printaspour(const gen & feuille,const string & sommetstr){
    if ( (feuille.type!=_VECT) || (feuille._VECTptr->size()!=4) )
      return sommetstr+'('+feuille.print()+')';
    const_iterateur it=feuille._VECTptr->begin(),itend=feuille._VECTptr->end();
    string res;
    if (is_zero(*it) && is_zero(*(it+2))){
      ++it;
      res ="tantque " + sametoequal(*it).print() + " faire ";
      ++it;
      ++it;
      indent_spaces += 2;
      if ((it->type==_SYMB) && (it->_SYMBptr->sommet==at_bloc))
          res += printasinnerbloc(it->_SYMBptr->feuille);
      else
          res += it->print() ;
      indent_spaces -= 2;
      return res+indent()+" ftantque;";
    }
    else {  
      if ( (it->type!=_SYMB) || (it->_SYMBptr->sommet!=at_sto) || ((it+2)->type!=_SYMB) || ((it+2)->_SYMBptr->sommet!=at_sto) || (it->_SYMBptr->feuille._VECTptr->back()!=(it+2)->_SYMBptr->feuille._VECTptr->back()) ){
	res="pour (";
	res += it->print() + ';';
	++it;
	res += it->print() + ';';
	++it;
	res += it->print() + ") ";
	++it;
	indent_spaces += 2;
	res += it->print() ;
	indent_spaces -= 2;
      }
      else {
	gen var_name=it->_SYMBptr->feuille._VECTptr->back();
	gen step=normal((it+2)->_SYMBptr->feuille._VECTptr->front()-var_name);
	gen condition=*(it+1),limite;
	bool simple_loop=false,strict=true,ascending=true;
	if (condition.type==_SYMB){
	  unary_function_ptr op=condition._SYMBptr->sommet;
	  if (condition._SYMBptr->feuille.type==_VECT){
	    if (op==at_inferieur)
	      simple_loop=true;
	    if (op==at_inf_ou_egal){
	      strict=false;
	      simple_loop=true;
	    }
	    if (op==at_superieur){
	      simple_loop=(maple_mode==2);
	      ascending=false;
	    }
	    if (op==at_sup_ou_egal){
	      simple_loop=(maple_mode==2);
	      ascending=false;
	      strict=false;
	    }
	  }
	  if (simple_loop){
	    simple_loop=(condition._SYMBptr->feuille._VECTptr->front()==var_name);
	    limite=condition._SYMBptr->feuille._VECTptr->back();
	  }
	}
	res ="pour ";
	res += var_name.print();
	res += " de ";
	res += it->_SYMBptr->feuille._VECTptr->front().print();
	if (simple_loop){
	  step = abs(step);
	  if (ascending)
	    res += " jusque ";
	  else
	    res += " downto ";
	  res += limite.print();
	  if (strict){
	    if (ascending)
	      res +="+";
	    else
	      res += "-";
	    res += step.print();
	    res += "/2";
	  }
	}
	if (!is_one(step)){
	  res += " by ";
	  res += step.print();
	}
	if (!simple_loop){
	  res += " tantque ";
	  res += (it+1)->print();
	}
	res += " faire ";
	it += 3;
	if ((it->type==_SYMB) && (it->_SYMBptr->sommet==at_bloc))
	  res += printasinnerbloc(it->_SYMBptr->feuille);
	else
	  res += it->print() ;
	return res + indent()+" fpour;";
      }
    }
    if (res[res.size()-1]!='}')
      res += "; ";
    return res;
  }
  gen _pour(const gen & g){
    return _for(g);
  }
  const string _pour_s("pour");
  unary_function_unary __pour(&_pour,_pour_s,&printaspour);
  unary_function_ptr at_pour (&__pour,_QUOTE_ARGUMENTS,T_FOR);

  const string _fpour_s("fpour");
  unary_function_unary __fpour(&_pour,_fpour_s);
  unary_function_ptr at_fpour (&__fpour,_QUOTE_ARGUMENTS,T_BLOC_END);
  
  const string _ftantque_s("ftantque");
  unary_function_unary __ftantque(&_pour,_ftantque_s);
  unary_function_ptr at_ftantque (&__ftantque,_QUOTE_ARGUMENTS,T_BLOC_END);
  
  const string _de_s("de");
  unary_function_unary __de(&_pour,_de_s);
  unary_function_ptr at_de (&__de,_QUOTE_ARGUMENTS,T_FROM);

  const string _faire_s("faire");
  unary_function_unary __faire(&_pour,_faire_s);
  unary_function_ptr at_faire (&__faire,_QUOTE_ARGUMENTS,T_DO);

  const string _pas_s("pas");
  unary_function_unary __pas(&_pour,_pas_s);
  unary_function_ptr at_pas (&__pas,_QUOTE_ARGUMENTS,T_BY);

  const string _jusque_s("jusque");
  unary_function_unary __jusque(&_pour,_jusque_s);
  unary_function_ptr at_jusque (&__jusque,_QUOTE_ARGUMENTS,T_TO);

  const string _tantque_s("tantque");
  unary_function_unary __tantque(&_pour,_tantque_s);
  unary_function_ptr at_tantque (&__tantque,_QUOTE_ARGUMENTS,T_WHILE);

  const string _et_s("et");
  unary_function_unary __et(&_and,_et_s);
  unary_function_ptr at_et (&__et,_QUOTE_ARGUMENTS,T_AND_OP);

  const string _oufr_s("ou");
  unary_function_unary __oufr(&_ou,_oufr_s);
  unary_function_ptr at_oufr (&__oufr,_QUOTE_ARGUMENTS,T_OR_OP);

  const string _non_s("non");
  unary_function_unary __non(&_not,_non_s);
  unary_function_ptr at_non (&__non,_QUOTE_ARGUMENTS,T_NOT);

  const string _resultat_s("resultat");
  unary_function_unary __resultat(&_nop,_resultat_s);
  unary_function_ptr at_resultat (&__resultat,_QUOTE_ARGUMENTS,T_RETURN);

#ifndef NO_NAMESPACE_GIAC
} // namespace giac
#endif // ndef NO_NAMESPACE_GIAC
