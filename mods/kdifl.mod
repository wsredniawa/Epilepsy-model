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
	  kbath   =  0.002 (mM)        : seawater (squid axon!)
	  fhspace = 50 (angstrom)  : effective thickness of F-H space
	  txfer   =  2 (ms)  : tau for F-H space <-> bath exchange 
	  difc = 3
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