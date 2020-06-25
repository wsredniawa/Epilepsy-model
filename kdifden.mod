    NEURON {
      SUFFIX kdifden
      USEION k READ ko WRITE ko
	  GLOBAL difsom,difax
	  POINTER kog
    }

    UNITS {
      (mV) = (millivolt)
      (mA) = (milliamp)
      FARADAY = (faraday) (coulombs)
      (molar) = (1/liter)
      (mM) = (millimolar)
    }

    PARAMETER {
	  difsom   =  1    : seawater (squid axon!)
	  difax  = 1
    }

    ASSIGNED {
	  kog
    }
    INITIAL {
	 ko = 3}
	 
	STATE {
     ko (mM)
	 :ki (mM)
    }

    BREAKPOINT {
       SOLVE state METHOD derivimplicit
	}
    DERIVATIVE state {
		ko'= -(ko-3)/difsom  + kog*difax
    }
