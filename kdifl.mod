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
	  crossdiff=1
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
		ko'=(1e6)*(ik)/(fhspace*FARADAY)-glia(difc, txfer, ko)-(ko-3)*kbath :+crossdiff*(koloc-ko)
    }
    FUNCTION glia(loc1, loc2, loc3) {
    glia = loc1/(1+pow(loc2,(-1*(loc3-15))))
    }