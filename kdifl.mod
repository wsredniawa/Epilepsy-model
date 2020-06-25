    NEURON {
      SUFFIX kdifl
      USEION k READ ik,ko WRITE ko:,k
      RANGE fhspace,txfer
	  GLOBAL kbath,difc
	  
    }

    UNITS {
      (mV) = (millivolt)
      (mA) = (milliamp)
      FARADAY = (faraday) (coulombs)
      (molar) = (1/liter)
      (mM) = (millimolar)
    }

    PARAMETER {
	  kbath   =  3 (mM)        : seawater (squid axon!)
	  fhspace = 300 (angstrom)  : effective thickness of F-H space
	  txfer   =  50 (ms)  : tau for F-H space <-> bath exchange 
	  difc = 0.0001
    }

    ASSIGNED {
      ik (mA/cm2)
    }

    STATE {
     ko (mM)
	 :ki (mM)
    }

    BREAKPOINT {
       SOLVE state METHOD derivimplicit
    }

    DERIVATIVE state {
		ko'=(1e6)*(ik)/(fhspace*FARADAY) - difc/(1+pow(txfer,(-1*(ko-15))))-(ko-3)*kbath
    }