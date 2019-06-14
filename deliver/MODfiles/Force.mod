TITLE Force.mod 

COMMENT
    tchuuuru
ENDCOMMENT

NEURON {
        POINT_PROCESS Force
        RANGE muforce, vpeak, v1, dv, tau
}

UNITS {
        (nA) = (nanoamp)
        (mV) = (millivolt)
       }

PARAMETER {
        del = 10 (ms)
        muforce = 0 (newton)
        dur = 100 (ms)
        tau = .1 (ms) <1e-3,1e6>
}

ASSIGNED {
        v (mV)
        i (nA)
        v1 (mV)
        vpeak (mV)
        dv (mV)
}

INITIAL {
         v1 = 0 (mV)
         vpeak = 0 (mV)
         dv = 0 (mV)
         }
         
BREAKPOINT {
        SOLVE check METHOD after_cvode
        at_time(del)
        at_time(del+dur)
        
        VERBATIM
        
        if(t<del){
            muforce = 0;
        }
        else{
             if(t<(del+dur)){
                muforce = v*2;
                muforce = vpeak * alpha( ((t-del)/tau) );
             }
             else{
                  muforce = 0;
             }
        }
            
        ENDVERBATIM
}

FUNCTION alpha(x) {
        if (x<0) {
            alpha = 0
        }else{
            alpha = x * exp(1-x)
        }
}
PROCEDURE check() {
        if (v > vpeak){
            vpeak = v
        }
        dv = v - v1
        v1 = v
}


