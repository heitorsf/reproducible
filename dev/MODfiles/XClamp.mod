TITLE XClamp.mod 

COMMENT
Since this is an electrode current, positive values of i depolarize the cell
and in the presence of the extracellular mechanism there will be a change
in vext since i is not a transmembrane current but a current injected
directly to the inside of the cell.

---------------- Signal types ---------------
Type #  Description;                Inputs
1       IClamp (NEURON original);   del, dur, amp
2       iTriang;                    del, predur, posdur, posamp
3       ISine;                      del, dur, amp, frq, phase, offset, sinedelay
4       IRand;                      del, dur, amp

ENDCOMMENT

NEURON {
        POINT_PROCESS IClamp_X
        RANGE   sig_type, del, dur, amp,
                predur, posdur, posamp,
                frq, phase, offset, sinedelay
        ELECTRODE_CURRENT i
}

UNITS {
        (nA) = (nanoamp)
             }

PARAMETER {
        sig_type = 1
    
        del = 10 (ms)
        dur = 50 (ms)
        amp = 1 (nA)   
        
        predur = 200 (ms)
        posdur = 200 (ms)
        posamp = 1 (nA)
        
        frq = 2 (Hz)
        phase = 0 (degree)
        offset = 0 (nA)
        sinedelay = 0 (ms)
}

ASSIGNED {
        i (nA)
}

BREAKPOINT {
        VERBATIM
        int typeint = sig_type;
        int phase_rad = 0;
        switch(typeint){

            case 1: // step
                at_time(nrn_threads, del);
                at_time(nrn_threads, del+dur);
                if (t < del + dur && t >= del)
                    i = amp;
                else
                    i = 0;
                break;                 
                 
            case 2: // triangular
                at_time(nrn_threads, del);
                at_time(nrn_threads, del + predur);
                at_time(nrn_threads, del + predur + posdur);
                if (t < del)
                    i = 0;
                else{
                    if (t < del+predur)
                        i = posamp * (t-del)/predur;
                    else{
                        if (t < del+predur+posdur)
                            i = posamp - (posamp * (t-(del+predur))/posdur);
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
                if(t < del) i = 0;
                else
                    if(t < del + dur) i = amp * rand()/RAND_MAX;
                    else
                        i = 0;
                break;
        }
        ENDVERBATIM
}
