#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;

extern void _Force_reg(void);
extern void _IClampX_reg(void);
extern void _XClamp_reg(void);
extern void _caL_reg(void);
extern void _force_syn_reg(void);
extern void _napp_reg(void);
extern void _peakcount_reg(void);
extern void _refrac_reg(void);

void modl_reg(){
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");

    fprintf(stderr," MODfiles//Force.mod");
    fprintf(stderr," MODfiles//IClampX.mod");
    fprintf(stderr," MODfiles//XClamp.mod");
    fprintf(stderr," MODfiles//caL.mod");
    fprintf(stderr," MODfiles//force_syn.mod");
    fprintf(stderr," MODfiles//napp.mod");
    fprintf(stderr," MODfiles//peakcount.mod");
    fprintf(stderr," MODfiles//refrac.mod");
    fprintf(stderr, "\n");
  }
  _Force_reg();
  _IClampX_reg();
  _XClamp_reg();
  _caL_reg();
  _force_syn_reg();
  _napp_reg();
  _peakcount_reg();
  _refrac_reg();
}
