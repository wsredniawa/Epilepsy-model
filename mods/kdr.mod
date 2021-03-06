COMMENT

	Maciej T. Lazarewicz, mlazarew@seas.upenn.edu

ENDCOMMENT


NEURON {

	SUFFIX kdr
	USEION k WRITE ik
	RANGE gkbar, ik
}
	
UNITS {

	(mA) = (milliamp)
	(mV) = (millivolt)
	(mS) = (millisiemens)
}

PARAMETER {

    gkbar =  15 (mS/cm2)
    ek   = -75 (mV)
}
    
ASSIGNED {

    v       (mV)
    ik      (mA/cm2)
    ninf    (1)
    taun    (ms)
}

STATE { n }

INITIAL { 

    rates(v)
    n  = ninf
}

BREAKPOINT {

	SOLVE states METHOD cnexp
	
	ik = (1e-3) * gkbar * n * (v-ek)
}


DERIVATIVE states { 

    rates(v)
    n' = (ninf-n)/taun
}

PROCEDURE rates(v(mV)) { LOCAL a, b

    a = fun3(v,  -24.9, -0.016,   -5)
    b = fun1(v,  -40,    0.25,   -40)
    
    ninf = a/(a+b)
    taun = 1.0/(a+b)
}

COMMENT

	Maciej T. Lazarewicz, mlazarew@seas.upenn.edu

ENDCOMMENT



:-------------------------------------------------------------------
FUNCTION fun1(v(mV),V0(mV),A(/ms),B(mV))(/ms) {

	 fun1 = A*exp((v-V0)/B)
}

FUNCTION fun2(v(mV),V0(mV),A(/ms),B(mV))(/ms) {

	 fun2 = A/(exp((v-V0)/B)+1)
}

FUNCTION fun3(v(mV),V0(mV),A(/ms),B(mV))(/ms) {

    if(fabs((v-V0)/B)<1e-6) {
    :if(v==V0) {
        fun3 = A*B/1(mV) * (1- 0.5 * (v-V0)/B)
    } else {
        fun3 = A/1(mV)*(v-V0)/(exp((v-V0)/B)-1)
    }
}

FUNCTION min(x,y) { if (x<=y){ min = x }else{ min = y } }





