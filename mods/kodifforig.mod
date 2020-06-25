    NEURON {
      SUFFIX kdiff2
      USEION k READ ik,ko WRITE ko:,ki
      RANGE fhspace,txfer
	  GLOBAL difsom, difax, difc
      :THREADSAFE
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
	  difsom   =  20    : seawater (squid axon!)
	  difax  = 20
	  fhspace = 300 (angstrom)  : effective thickness of F-H space
	  txfer   =  50 (ms)  : tau for F-H space <-> bath exchange 
	  difc = 0.0001
    }

    ASSIGNED {
      ik (mA/cm2)
      kog
    }
    
	STATE {
     ko (mM)
	 :ki (mM)
    }

    BREAKPOINT {
       SOLVE state METHOD derivimplicit
	}
    DERIVATIVE state {
		ko'=(1e6)*(ik)/(fhspace*FARADAY) - difc/(1+pow(txfer,(-1*(ko-15))))-(ko-3)*difsom+(kog - ko)*difax
    }
:ustawic sigmoide na 15 k2 = 0.0008/(1+exp((ko-15)/-1.09)

