/* Created by Language version: 6.2.0 */
/* NOT VECTORIZED */
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
 
#define _threadargscomma_ /**/
#define _threadargs_ /**/
 
#define _threadargsprotocomma_ /**/
#define _threadargsproto_ /**/
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 static double *_p; static Datum *_ppvar;
 
#define t nrn_threads->_t
#define dt nrn_threads->_dt
#define signal_type _p[0]
#define del _p[1]
#define dur _p[2]
#define amp _p[3]
#define risedur _p[4]
#define falldur _p[5]
#define after_fall_amp _p[6]
#define frq _p[7]
#define phase _p[8]
#define offset _p[9]
#define sinedelay _p[10]
#define mu _p[11]
#define sigma _p[12]
#define std _p[13]
#define media _p[14]
#define i _p[15]
#define ibb _p[16]
#define _g _p[17]
#define _nd_area  *_ppvar[0]._pval
 
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
 /* external NEURON variables */
 /* declaration of user functions */
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 extern Prop* nrn_point_prop_;
 static int _pointtype;
 static void* _hoc_create_pnt(_ho) Object* _ho; { void* create_point_process();
 return create_point_process(_pointtype, _ho);
}
 static void _hoc_destroy_pnt();
 static double _hoc_loc_pnt(_vptr) void* _vptr; {double loc_point_process();
 return loc_point_process(_pointtype, _vptr);
}
 static double _hoc_has_loc(_vptr) void* _vptr; {double has_loc_point();
 return has_loc_point(_vptr);
}
 static double _hoc_get_loc_pnt(_vptr)void* _vptr; {
 double get_loc_point_process(); return (get_loc_point_process(_vptr));
}
 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _p = _prop->param; _ppvar = _prop->dparam;
 }
 static void _hoc_setdata(void* _vptr) { Prop* _prop;
 _prop = ((Point_process*)_vptr)->_prop;
   _setdata(_prop);
 }
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 0,0
};
 static Member_func _member_func[] = {
 "loc", _hoc_loc_pnt,
 "has_loc", _hoc_has_loc,
 "get_loc", _hoc_get_loc_pnt,
 0, 0
};
 /* declare global and static user variables */
#define high high_IClampX
 double high = 1;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "del", "ms",
 "dur", "ms",
 "amp", "nA",
 "risedur", "ms",
 "falldur", "ms",
 "after_fall_amp", "nA",
 "frq", "Hz",
 "phase", "degree",
 "offset", "nA",
 "sinedelay", "ms",
 "mu", "nA",
 "sigma", "nA",
 "i", "nA",
 0,0
};
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "high_IClampX", &high_IClampX,
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
 static void _hoc_destroy_pnt(_vptr) void* _vptr; {
   destroy_point_process(_vptr);
}
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "6.2.0",
"IClampX",
 "signal_type",
 "del",
 "dur",
 "amp",
 "risedur",
 "falldur",
 "after_fall_amp",
 "frq",
 "phase",
 "offset",
 "sinedelay",
 "mu",
 "sigma",
 "std",
 "media",
 0,
 "i",
 0,
 0,
 0};
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
  if (nrn_point_prop_) {
	_prop->_alloc_seq = nrn_point_prop_->_alloc_seq;
	_p = nrn_point_prop_->param;
	_ppvar = nrn_point_prop_->dparam;
 }else{
 	_p = nrn_prop_data_alloc(_mechtype, 18, _prop);
 	/*initialize range parameters*/
 	signal_type = 1;
 	del = 10;
 	dur = 50;
 	amp = 1;
 	risedur = 200;
 	falldur = 200;
 	after_fall_amp = 0;
 	frq = 2;
 	phase = 0;
 	offset = 0;
 	sinedelay = 0;
 	mu = 0;
 	sigma = 1;
 	std = 1;
 	media = 0;
  }
 	_prop->param = _p;
 	_prop->param_size = 18;
  if (!nrn_point_prop_) {
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 2, _prop);
  }
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 
}
 static void _initlists();
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*f)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _IClampX_reg() {
	int _vectorized = 0;
  _initlists();
 	_pointtype = point_register_mech(_mechanism,
	 nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init,
	 hoc_nrnpointerindex, 0,
	 _hoc_create_pnt, _hoc_destroy_pnt, _member_func);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
  hoc_register_prop_size(_mechtype, 18, 2);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 IClampX /home/heitor/Dropbox/Unicamp/reproducible/dev/x86_64/IClampX.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "IClampX.mod ";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
 
/*VERBATIM*/
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <unistd.h>

static void initmodel() {
  int _i; double _save;_ninits++;
{
 {
   
/*VERBATIM*/
        #include <time.h>
        srand(getpid()*time(NULL));
 }

}
}

static void nrn_init(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if EXTRACELLULAR
 _nd = _ml->_nodelist[_iml];
 if (_nd->_extnode) {
    _v = NODEV(_nd) +_nd->_extnode->_v[0];
 }else
#endif
 {
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 }
 v = _v;
 initmodel();
}}

static double _nrn_current(double _v){double _current=0.;v=_v;{ {
   
/*VERBATIM*/
                
        int typeint = signal_type;
        int phase_rad = 0;
        switch(typeint){

            case 1: // step, NEURON's IClamp
                at_time(nrn_threads, del);
                at_time(nrn_threads, del+dur);
                if (t < del + dur && t >= del)
                    i = amp;
                else
                    i = 0;
                break;                 
                 
            case 2: // triangular
                at_time(nrn_threads, del);
                at_time(nrn_threads, del + risedur);
                at_time(nrn_threads, del + risedur + falldur);
                if (t < del)
                    i = 0;
                else{
                    if (t < del+risedur)
                        i = amp * (t-del)/risedur;
                    else{
                        if (t < del+risedur+falldur)
                            i = amp - ((amp-after_fall_amp) * (t-(del+risedur))/falldur);
                        else
                            i = 0;
                    }
                }
                break;
                
            case 3: // sine
                at_time(nrn_threads, del);
                at_time(nrn_threads, del + sinedelay);
                at_time(nrn_threads, del + dur);
                //phase_rad = M_PI*phase/180;
                if(t < del)
                    i = 0;
                else{
                    if(t < del + sinedelay)
                        i = offset;
                    else{
                        if(t < del + dur)
                            i = offset + amp*sin(2*M_PI*(frq/1000)*(t-del-sinedelay)
                                                     + (M_PI*phase/180));
                        else
                            i = 0;
                    }
                }
                break;
            
            case 4: // random
                at_time(nrn_threads, del);
                at_time(nrn_threads, del + dur);
                double a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12;
                double sum;
                if(t < del){
                     i = 0;
                }
                else{
                    if(t < del + dur){
                        //a1 =  rand()/randmax;
                        //a2 =  rand()/randmax;
                        //a3 =  rand()/randmax;
                        //a4 =  rand()/randmax;
                        //a5 =  rand()/randmax;
                        //a6 =  rand()/randmax;
                        //a7 =  rand()/randmax;
                        //a8 =  rand()/randmax;
                        //a9 =  rand()/randmax;
                        //a10 = rand()/randmax;
                        //a11 = rand()/randmax;
                        //a12 = rand()/randmax;
                        //sum = a1+a2+a3+a4+a5+a6+a7+a8+a9+a10+a11+a12;
                        //i = std*(sum-6) + media;

                        i = ibb;
                        
                    }
                    else
                        i = 0;
                }
                break;
        }
 }
 _current += i;

} return _current;
}

static void nrn_cur(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if EXTRACELLULAR
 _nd = _ml->_nodelist[_iml];
 if (_nd->_extnode) {
    _v = NODEV(_nd) +_nd->_extnode->_v[0];
 }else
#endif
 {
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 }
 _g = _nrn_current(_v + .001);
 	{ _rhs = _nrn_current(_v);
 	}
 _g = (_g - _rhs)/.001;
 _g *=  1.e2/(_nd_area);
 _rhs *= 1.e2/(_nd_area);
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) += _rhs;
  }else
#endif
  {
	NODERHS(_nd) += _rhs;
  }
  if (_nt->_nrn_fast_imem) { _nt->_nrn_fast_imem->_nrn_sav_rhs[_ni[_iml]] += _rhs; }
#if EXTRACELLULAR
 if (_nd->_extnode) {
   *_nd->_extnode->_rhs[0] += _rhs;
 }
#endif
 
}}

static void nrn_jacob(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) -= _g;
  }else
#endif
  {
	NODED(_nd) -= _g;
  }
  if (_nt->_nrn_fast_imem) { _nt->_nrn_fast_imem->_nrn_sav_d[_ni[_iml]] -= _g; }
#if EXTRACELLULAR
 if (_nd->_extnode) {
   *_nd->_extnode->_d[0] += _g;
 }
#endif
 
}}

static void nrn_state(_NrnThread* _nt, _Memb_list* _ml, int _type){

}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
_first = 0;
}
