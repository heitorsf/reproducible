NEURON {
	POINT_PROCESS PeakCount
	RANGE n, thresh, time, firing
	THREADSAFE : if APCount.record uses distinct instances of Vector
}

UNITS {
	(mV) = (millivolt)
}

PARAMETER {
	n
	thresh = 60 (mV)
	time (ms)
    vm1 = 0
    vm2 = 0
    tm1 = 0
}

ASSIGNED {
	firing
	space
}

VERBATIM
extern void vector_resize();
extern double* vector_vec();
extern void* vector_arg();
ENDVERBATIM

INITIAL {
	n = 0
	firing = 0
VERBATIM
	{ void* vv;
		vv = *((void**)(&space));
		if (vv) {
			vector_resize(vv, 0);
		}
	}
ENDVERBATIM
    net_send(firing,3)
	check()
}

BREAKPOINT {
	SOLVE check METHOD after_cvode
}

PROCEDURE check() {
VERBATIM
    if(vm1 > v && vm1 > vm2 && vm1 > thresh && !firing){
		int size; double* px; void* vv;
		firing = 1;
		time = tm1;
		n += 1.;
		vv = *((void**)(&space));
		if (vv) {
			size = (int)n;
			vector_resize(vv, size);
			px = vector_vec(vv);
			px[size-1] = time;
		}
	}
	if (firing && v < thresh && t > time) {
		firing = 0;
	}
    vm2 = vm1;
    vm1 = v;
    tm1 = t;
ENDVERBATIM
}

PROCEDURE record() {
VERBATIM
	extern void* vector_arg();
	void** vv;
	vv = (void**)(&space);
	*vv = (void*)0;
	if (ifarg(1)) {
		*vv = vector_arg(1);
	}
ENDVERBATIM
}

NET_RECEIVE (w) {
    if (flag == 3){
        WATCH (firing > 0) 2
    }
    if (flag == 2){
        net_event(t)
    }
}

