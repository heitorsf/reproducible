/* Created by Language version: 6.2.0 */
/* VECTORIZED */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "scoplib_ansi.h"
#undef PI
#define nil 0
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define _threadargscomma_ _p, _ppvar, _thread, _nt,
#define _threadargs_ _p, _ppvar, _thread, _nt
 
#define _threadargsprotocomma_ double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt,
#define _threadargsproto_ double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 /* Thread safe. No static _p or _ppvar. */
 
#define t _nt->_t
#define dt _nt->_dt
#define gnabar _p[0]
#define gnapbar _p[1]
#define gkfbar _p[2]
#define gksbar _p[3]
#define gl _p[4]
#define el _p[5]
#define vtraub _p[6]
#define mact _p[7]
#define rinact _p[8]
#define gna _p[9]
#define gnap _p[10]
#define gkf _p[11]
#define gks _p[12]
#define inap _p[13]
#define inaf _p[14]
#define ikf _p[15]
#define iks _p[16]
#define il _p[17]
#define m _p[18]
#define h _p[19]
#define p _p[20]
#define n _p[21]
#define r _p[22]
#define Dm _p[23]
#define Dh _p[24]
#define Dp _p[25]
#define Dn _p[26]
#define Dr _p[27]
#define ena _p[28]
#define ek _p[29]
#define ina _p[30]
#define ik _p[31]
#define v _p[32]
#define _g _p[33]
#define _ion_ena	*_ppvar[0]._pval
#define _ion_ina	*_ppvar[1]._pval
#define _ion_dinadv	*_ppvar[2]._pval
#define _ion_ek	*_ppvar[3]._pval
#define _ion_ik	*_ppvar[4]._pval
#define _ion_dikdv	*_ppvar[5]._pval
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  -1;
 static Datum* _extcall_thread;
 static Prop* _extcall_prop;
 /* external NEURON variables */
 extern double celsius;
 /* declaration of user functions */
 static void _hoc_rates(void);
 static void _hoc_vtrap(void);
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _extcall_prop = _prop;
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 "setdata_napp", _hoc_setdata,
 "rates_napp", _hoc_rates,
 "vtrap_napp", _hoc_vtrap,
 0, 0
};
#define vtrap vtrap_napp
 extern double vtrap( _threadargsprotocomma_ double , double );
 
static void _check_rates(double*, Datum*, Datum*, _NrnThread*); 
static void _check_table_thread(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, int _type) {
   _check_rates(_p, _ppvar, _thread, _nt);
 }
 /* declare global and static user variables */
 static int _thread1data_inuse = 0;
static double _thread1data[10];
#define _gth 0
#define htau_napp _thread1data[0]
#define htau _thread[_gth]._pval[0]
#define hinf_napp _thread1data[1]
#define hinf _thread[_gth]._pval[1]
#define mtau_napp _thread1data[2]
#define mtau _thread[_gth]._pval[2]
#define minf_napp _thread1data[3]
#define minf _thread[_gth]._pval[3]
#define ntau_napp _thread1data[4]
#define ntau _thread[_gth]._pval[4]
#define ninf_napp _thread1data[5]
#define ninf _thread[_gth]._pval[5]
#define ptau_napp _thread1data[6]
#define ptau _thread[_gth]._pval[6]
#define pinf_napp _thread1data[7]
#define pinf _thread[_gth]._pval[7]
#define rtau_napp _thread1data[8]
#define rtau _thread[_gth]._pval[8]
#define rinf_napp _thread1data[9]
#define rinf _thread[_gth]._pval[9]
#define usetable usetable_napp
 double usetable = 1;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "gl_napp", 0, 1e+09,
 "gksbar_napp", 0, 1e+09,
 "gkfbar_napp", 0, 1e+09,
 "gnapbar_napp", 0, 1e+09,
 "gnabar_napp", 0, 1e+09,
 "usetable_napp", 0, 1,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "mtau_napp", "ms",
 "htau_napp", "ms",
 "ptau_napp", "ms",
 "ntau_napp", "ms",
 "rtau_napp", "ms",
 "gnabar_napp", "S/cm2",
 "gnapbar_napp", "S/cm2",
 "gkfbar_napp", "S/cm2",
 "gksbar_napp", "S/cm2",
 "gl_napp", "S/cm2",
 "el_napp", "mV",
 "vtraub_napp", "mV",
 "mact_napp", "mV",
 "rinact_napp", "/ms",
 "gna_napp", "S/cm2",
 "gnap_napp", "S/cm2",
 "gkf_napp", "S/cm2",
 "gks_napp", "S/cm2",
 "inap_napp", "mA/cm2",
 "inaf_napp", "mA/cm2",
 "ikf_napp", "mA/cm2",
 "iks_napp", "mA/cm2",
 "il_napp", "mA/cm2",
 0,0
};
 static double delta_t = 0.01;
 static double h0 = 0;
 static double m0 = 0;
 static double n0 = 0;
 static double p0 = 0;
 static double r0 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "minf_napp", &minf_napp,
 "hinf_napp", &hinf_napp,
 "pinf_napp", &pinf_napp,
 "ninf_napp", &ninf_napp,
 "rinf_napp", &rinf_napp,
 "mtau_napp", &mtau_napp,
 "htau_napp", &htau_napp,
 "ptau_napp", &ptau_napp,
 "ntau_napp", &ntau_napp,
 "rtau_napp", &rtau_napp,
 "usetable_napp", &usetable_napp,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(_NrnThread*, _Memb_list*, int);
static void nrn_state(_NrnThread*, _Memb_list*, int);
 static void nrn_cur(_NrnThread*, _Memb_list*, int);
static void  nrn_jacob(_NrnThread*, _Memb_list*, int);
 
static int _ode_count(int);
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(_NrnThread*, _Memb_list*, int);
static void _ode_matsol(_NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[6]._i
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "6.2.0",
"napp",
 "gnabar_napp",
 "gnapbar_napp",
 "gkfbar_napp",
 "gksbar_napp",
 "gl_napp",
 "el_napp",
 "vtraub_napp",
 "mact_napp",
 "rinact_napp",
 0,
 "gna_napp",
 "gnap_napp",
 "gkf_napp",
 "gks_napp",
 "inap_napp",
 "inaf_napp",
 "ikf_napp",
 "iks_napp",
 "il_napp",
 0,
 "m_napp",
 "h_napp",
 "p_napp",
 "n_napp",
 "r_napp",
 0,
 0};
 static Symbol* _na_sym;
 static Symbol* _k_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 34, _prop);
 	/*initialize range parameters*/
 	gnabar = 0.03;
 	gnapbar = 3.3e-05;
 	gkfbar = 0.016;
 	gksbar = 0.004;
 	gl = 0.0003;
 	el = -54.3;
 	vtraub = 50;
 	mact = 15;
 	rinact = 0.05;
 	_prop->param = _p;
 	_prop->param_size = 34;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 7, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_na_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[0]._pval = &prop_ion->param[0]; /* ena */
 	_ppvar[1]._pval = &prop_ion->param[3]; /* ina */
 	_ppvar[2]._pval = &prop_ion->param[4]; /* _ion_dinadv */
 prop_ion = need_memb(_k_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[3]._pval = &prop_ion->param[0]; /* ek */
 	_ppvar[4]._pval = &prop_ion->param[3]; /* ik */
 	_ppvar[5]._pval = &prop_ion->param[4]; /* _ion_dikdv */
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 static void _thread_mem_init(Datum*);
 static void _thread_cleanup(Datum*);
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*f)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _napp_reg() {
	int _vectorized = 1;
  _initlists();
 	ion_reg("na", -10000.);
 	ion_reg("k", -10000.);
 	_na_sym = hoc_lookup("na_ion");
 	_k_sym = hoc_lookup("k_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 2);
  _extcall_thread = (Datum*)ecalloc(1, sizeof(Datum));
  _thread_mem_init(_extcall_thread);
  _thread1data_inuse = 0;
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 1, _thread_mem_init);
     _nrn_thread_reg(_mechtype, 0, _thread_cleanup);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
     _nrn_thread_table_reg(_mechtype, _check_table_thread);
  hoc_register_prop_size(_mechtype, 34, 7);
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 napp /home/heitor/Dropbox/Unicamp/reproducible/dev/x86_64/napp.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double *_t_minf;
 static double *_t_mtau;
 static double *_t_hinf;
 static double *_t_htau;
 static double *_t_pinf;
 static double *_t_ptau;
 static double *_t_ninf;
 static double *_t_ntau;
 static double *_t_rinf;
 static double *_t_rtau;
static int _reset;
static char *modelname = "napp.mod   squid sodium, potassium, and leak channels";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int _f_rates(_threadargsprotocomma_ double);
static int rates(_threadargsprotocomma_ double);
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static void _n_rates(_threadargsprotocomma_ double _lv);
 static int _slist1[5], _dlist1[5];
 static int states(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {int _reset = 0; {
   rates ( _threadargscomma_ v ) ;
   Dm = ( minf - m ) / mtau ;
   Dh = ( hinf - h ) / htau ;
   Dp = ( pinf - p ) / ptau ;
   Dn = ( ninf - n ) / ntau ;
   Dr = ( rinf - r ) / rtau ;
   }
 return _reset;
}
 static int _ode_matsol1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
 rates ( _threadargscomma_ v ) ;
 Dm = Dm  / (1. - dt*( ( ( ( - 1.0 ) ) ) / mtau )) ;
 Dh = Dh  / (1. - dt*( ( ( ( - 1.0 ) ) ) / htau )) ;
 Dp = Dp  / (1. - dt*( ( ( ( - 1.0 ) ) ) / ptau )) ;
 Dn = Dn  / (1. - dt*( ( ( ( - 1.0 ) ) ) / ntau )) ;
 Dr = Dr  / (1. - dt*( ( ( ( - 1.0 ) ) ) / rtau )) ;
 return 0;
}
 /*END CVODE*/
 static int states (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) { {
   rates ( _threadargscomma_ v ) ;
    m = m + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / mtau)))*(- ( ( ( minf ) ) / mtau ) / ( ( ( ( - 1.0) ) ) / mtau ) - m) ;
    h = h + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / htau)))*(- ( ( ( hinf ) ) / htau ) / ( ( ( ( - 1.0) ) ) / htau ) - h) ;
    p = p + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / ptau)))*(- ( ( ( pinf ) ) / ptau ) / ( ( ( ( - 1.0) ) ) / ptau ) - p) ;
    n = n + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / ntau)))*(- ( ( ( ninf ) ) / ntau ) / ( ( ( ( - 1.0) ) ) / ntau ) - n) ;
    r = r + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / rtau)))*(- ( ( ( rinf ) ) / rtau ) / ( ( ( ( - 1.0) ) ) / rtau ) - r) ;
   }
  return 0;
}
 static double _mfac_rates, _tmin_rates;
  static void _check_rates(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
  static int _maktable=1; int _i, _j, _ix = 0;
  double _xi, _tmax;
  static double _sav_celsius;
  if (!usetable) {return;}
  if (_sav_celsius != celsius) { _maktable = 1;}
  if (_maktable) { double _x, _dx; _maktable=0;
   _tmin_rates =  - 100.0 ;
   _tmax =  100.0 ;
   _dx = (_tmax - _tmin_rates)/200.; _mfac_rates = 1./_dx;
   for (_i=0, _x=_tmin_rates; _i < 201; _x += _dx, _i++) {
    _f_rates(_p, _ppvar, _thread, _nt, _x);
    _t_minf[_i] = minf;
    _t_mtau[_i] = mtau;
    _t_hinf[_i] = hinf;
    _t_htau[_i] = htau;
    _t_pinf[_i] = pinf;
    _t_ptau[_i] = ptau;
    _t_ninf[_i] = ninf;
    _t_ntau[_i] = ntau;
    _t_rinf[_i] = rinf;
    _t_rtau[_i] = rtau;
   }
   _sav_celsius = celsius;
  }
 }

 static int rates(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, double _lv) { 
#if 0
_check_rates(_p, _ppvar, _thread, _nt);
#endif
 _n_rates(_p, _ppvar, _thread, _nt, _lv);
 return 0;
 }

 static void _n_rates(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, double _lv){ int _i, _j;
 double _xi, _theta;
 if (!usetable) {
 _f_rates(_p, _ppvar, _thread, _nt, _lv); return; 
}
 _xi = _mfac_rates * (_lv - _tmin_rates);
 if (isnan(_xi)) {
  minf = _xi;
  mtau = _xi;
  hinf = _xi;
  htau = _xi;
  pinf = _xi;
  ptau = _xi;
  ninf = _xi;
  ntau = _xi;
  rinf = _xi;
  rtau = _xi;
  return;
 }
 if (_xi <= 0.) {
 minf = _t_minf[0];
 mtau = _t_mtau[0];
 hinf = _t_hinf[0];
 htau = _t_htau[0];
 pinf = _t_pinf[0];
 ptau = _t_ptau[0];
 ninf = _t_ninf[0];
 ntau = _t_ntau[0];
 rinf = _t_rinf[0];
 rtau = _t_rtau[0];
 return; }
 if (_xi >= 200.) {
 minf = _t_minf[200];
 mtau = _t_mtau[200];
 hinf = _t_hinf[200];
 htau = _t_htau[200];
 pinf = _t_pinf[200];
 ptau = _t_ptau[200];
 ninf = _t_ninf[200];
 ntau = _t_ntau[200];
 rinf = _t_rinf[200];
 rtau = _t_rtau[200];
 return; }
 _i = (int) _xi;
 _theta = _xi - (double)_i;
 minf = _t_minf[_i] + _theta*(_t_minf[_i+1] - _t_minf[_i]);
 mtau = _t_mtau[_i] + _theta*(_t_mtau[_i+1] - _t_mtau[_i]);
 hinf = _t_hinf[_i] + _theta*(_t_hinf[_i+1] - _t_hinf[_i]);
 htau = _t_htau[_i] + _theta*(_t_htau[_i+1] - _t_htau[_i]);
 pinf = _t_pinf[_i] + _theta*(_t_pinf[_i+1] - _t_pinf[_i]);
 ptau = _t_ptau[_i] + _theta*(_t_ptau[_i+1] - _t_ptau[_i]);
 ninf = _t_ninf[_i] + _theta*(_t_ninf[_i+1] - _t_ninf[_i]);
 ntau = _t_ntau[_i] + _theta*(_t_ntau[_i+1] - _t_ntau[_i]);
 rinf = _t_rinf[_i] + _theta*(_t_rinf[_i+1] - _t_rinf[_i]);
 rtau = _t_rtau[_i] + _theta*(_t_rtau[_i+1] - _t_rtau[_i]);
 }

 
static int  _f_rates ( _threadargsprotocomma_ double _lv ) {
   double _lalpha , _lbeta , _lsum , _lv2 ;
  _lv2 = _lv - vtraub ;
   _lalpha = .64 * vtrap ( _threadargscomma_ mact - _lv2 , 4.0 ) ;
   _lbeta = .56 * vtrap ( _threadargscomma_ _lv2 - 40.0 , 5.0 ) ;
   _lsum = _lalpha + _lbeta ;
   mtau = 1.0 / _lsum ;
   minf = _lalpha / _lsum ;
   _lalpha = .928 * exp ( ( 17.0 - _lv2 ) / 18.0 ) ;
   _lbeta = 9.0 / ( exp ( ( 40.0 - _lv2 ) / 5.0 ) + 1.0 ) ;
   _lsum = _lalpha + _lbeta ;
   htau = 1.0 / _lsum ;
   hinf = _lalpha / _lsum ;
   _lalpha = .64 * vtrap ( _threadargscomma_ 5.0 - _lv2 , 4.0 ) ;
   _lbeta = .56 * vtrap ( _threadargscomma_ _lv2 - 30.0 , 5.0 ) ;
   _lsum = _lalpha + _lbeta ;
   ptau = 1.0 / _lsum ;
   pinf = _lalpha / _lsum ;
   _lalpha = .08 * vtrap ( _threadargscomma_ 15.0 - _lv2 , 7.0 ) ;
   _lbeta = 2.0 * exp ( ( 10.0 - _lv2 ) / 40.0 ) ;
   _lsum = _lalpha + _lbeta ;
   ntau = 1.0 / _lsum ;
   ninf = _lalpha / _lsum ;
   _lalpha = ( 3.5 ) / ( exp ( ( 55.0 - _lv2 ) / 4.0 ) + 1.0 ) ;
   _lbeta = rinact ;
   _lsum = _lalpha + _lbeta ;
   rtau = 1.0 / _lsum ;
   rinf = _lalpha / _lsum ;
    return 0; }
 
static void _hoc_rates(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 
#if 1
 _check_rates(_p, _ppvar, _thread, _nt);
#endif
 _r = 1.;
 rates ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
double vtrap ( _threadargsprotocomma_ double _lx , double _ly ) {
   double _lvtrap;
 if ( fabs ( _lx / _ly ) < 1e-6 ) {
     _lvtrap = _ly * ( 1.0 - _lx / _ly / 2.0 ) ;
     }
   else {
     _lvtrap = _lx / ( exp ( _lx / _ly ) - 1.0 ) ;
     }
   
return _lvtrap;
 }
 
static void _hoc_vtrap(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r =  vtrap ( _p, _ppvar, _thread, _nt, *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ return 5;}
 
static void _ode_spec(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ena = _ion_ena;
  ek = _ion_ek;
     _ode_spec1 (_p, _ppvar, _thread, _nt);
   }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
	double* _p; Datum* _ppvar;
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 5; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ena = _ion_ena;
  ek = _ion_ek;
 _ode_matsol1 (_p, _ppvar, _thread, _nt);
 }}
 
static void _thread_mem_init(Datum* _thread) {
  if (_thread1data_inuse) {_thread[_gth]._pval = (double*)ecalloc(10, sizeof(double));
 }else{
 _thread[_gth]._pval = _thread1data; _thread1data_inuse = 1;
 }
 }
 
static void _thread_cleanup(Datum* _thread) {
  if (_thread[_gth]._pval == _thread1data) {
   _thread1data_inuse = 0;
  }else{
   free((void*)_thread[_gth]._pval);
  }
 }
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_na_sym, _ppvar, 0, 0);
   nrn_update_ion_pointer(_na_sym, _ppvar, 1, 3);
   nrn_update_ion_pointer(_na_sym, _ppvar, 2, 4);
   nrn_update_ion_pointer(_k_sym, _ppvar, 3, 0);
   nrn_update_ion_pointer(_k_sym, _ppvar, 4, 3);
   nrn_update_ion_pointer(_k_sym, _ppvar, 5, 4);
 }

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
  int _i; double _save;{
  h = h0;
  m = m0;
  n = n0;
  p = p0;
  r = r0;
 {
   rates ( _threadargscomma_ v ) ;
   m = minf ;
   h = hinf ;
   p = pinf ;
   n = ninf ;
   r = rinf ;
   }
 
}
}

static void nrn_init(_NrnThread* _nt, _Memb_list* _ml, int _type){
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];

#if 0
 _check_rates(_p, _ppvar, _thread, _nt);
#endif
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
  ena = _ion_ena;
  ek = _ion_ek;
 initmodel(_p, _ppvar, _thread, _nt);
  }}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   gna = gnabar * m * m * m * h ;
   gnap = gnapbar * p * p * p ;
   inaf = gna * ( v - ena ) ;
   inap = gnap * ( v - ena ) ;
   ina = inaf + inap ;
   gkf = gkfbar * n * n * n * n ;
   gks = gksbar * r * r ;
   ikf = ( gkf ) * ( v - ek ) ;
   iks = ( gks ) * ( v - ek ) ;
   ik = ikf + iks ;
   il = gl * ( v - el ) ;
   }
 _current += ina;
 _current += ik;
 _current += il;

} return _current;
}

static void nrn_cur(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
  ena = _ion_ena;
  ek = _ion_ek;
 _g = _nrn_current(_p, _ppvar, _thread, _nt, _v + .001);
 	{ double _dik;
 double _dina;
  _dina = ina;
  _dik = ik;
 _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
  _ion_dinadv += (_dina - ina)/.001 ;
  _ion_dikdv += (_dik - ik)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ina += ina ;
  _ion_ik += ik ;
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}}

static void nrn_jacob(_NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}}

static void nrn_state(_NrnThread* _nt, _Memb_list* _ml, int _type) {
 double _break, _save;
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 _break = t + .5*dt; _save = t;
 v=_v;
{
  ena = _ion_ena;
  ek = _ion_ek;
 { {
 for (; t < _break; t += dt) {
   states(_p, _ppvar, _thread, _nt);
  
}}
 t = _save;
 }  }}

}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(m) - _p;  _dlist1[0] = &(Dm) - _p;
 _slist1[1] = &(h) - _p;  _dlist1[1] = &(Dh) - _p;
 _slist1[2] = &(p) - _p;  _dlist1[2] = &(Dp) - _p;
 _slist1[3] = &(n) - _p;  _dlist1[3] = &(Dn) - _p;
 _slist1[4] = &(r) - _p;  _dlist1[4] = &(Dr) - _p;
   _t_minf = makevector(201*sizeof(double));
   _t_mtau = makevector(201*sizeof(double));
   _t_hinf = makevector(201*sizeof(double));
   _t_htau = makevector(201*sizeof(double));
   _t_pinf = makevector(201*sizeof(double));
   _t_ptau = makevector(201*sizeof(double));
   _t_ninf = makevector(201*sizeof(double));
   _t_ntau = makevector(201*sizeof(double));
   _t_rinf = makevector(201*sizeof(double));
   _t_rtau = makevector(201*sizeof(double));
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif
