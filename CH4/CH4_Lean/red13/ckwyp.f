C
C     A 13-species reduced mechanism for lean atmospheric CH4/air
C
C     Contact: 
C     Tianfeng Lu
C     University of Connecticut
C     Email: tlu@engr.uconn.edu
C
C     Reference:
C     R. Sankaran, E.R. Hawkes, J.H. Chen, T.F. Lu, C.K. Law, 
C     "Structure of a spatially developing turbulent lean methane–air Bunsen flame," 
C     Proceedings of the Combustion Institute 31 (2007) 1291–1298.
C
      SUBROUTINE CKWYP  (P, T, Y, ICKWRK, RCKWRK, WDOT)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
C
      PARAMETER (RU=8.314510D7, RUC=RU/4.184D7, PATM=1.01325D6)
      DIMENSION Y(*), ICKWRK(*), RCKWRK(*), WDOT(*)
      DIMENSION C(13), XQ(4) 
      DIMENSION RF(73), RB(73), ROP(73)
      DIMENSION EG(17), EQK(73), SMH(17), CTB(73), RKLOW(5)
      DATA RF, RB /73*0.D0, 73*0.D0/
C
	SMALL = 1.D-50
C
C	 convert Y to C
C
      C(1) = MAX(Y(1), SMALL)*4.96046521D-1
      C(2) = MAX(Y(2), SMALL)*9.92093043D-1
      C(3) = MAX(Y(3), SMALL)*6.25023433D-2
      C(4) = MAX(Y(4), SMALL)*3.12511716D-2
      C(5) = MAX(Y(5), SMALL)*5.87980383D-2
      C(6) = MAX(Y(6), SMALL)*5.55082499D-2
      C(7) = MAX(Y(7), SMALL)*3.02968146D-2
      C(8) = MAX(Y(8), SMALL)*6.65112065D-2
      C(9) = MAX(Y(9), SMALL)*6.23323639D-2
      C(10) = MAX(Y(10), SMALL)*3.57008335D-2
      C(11) = MAX(Y(11), SMALL)*2.27221341D-2
      C(12) = MAX(Y(12), SMALL)*3.33039255D-2
      C(13) = MAX(Y(13), SMALL)*3.56972032D-2
C
      SUM = 0.0
      DO K = 1, 13
         SUM = SUM + C(K)
      ENDDO
      SUM = P/(SUM*T*8.314510D7)
C
      DO K = 1, 13
         C(K) = C(K) * SUM
      ENDDO
C
C	 forward reaction rates
C
      ALOGT = LOG(T)
      TI = 1.0D0/T
      TI2 = TI*TI
C
      RF(1) = 1.2D17*TI
      RF(2) = 5.D17*TI
      RF(3) = EXP(1.08197783D1 +2.67D0*ALOGT -3.16523284D3*TI)
      RF(4) = 2.D13
      RF(5) = 8.D13
      RF(6) = 1.5D13
      RF(7) = 8.43D13
      RF(8) = EXP(2.07430685D1 +1.5D0*ALOGT -4.32766334D3*TI)
      RF(9) = EXP(3.40312786D1 -1.50965D3*TI)
      RF(10) = 3.D13
      RF(11) = 3.D13
      RF(12) = EXP(3.12945828D1 -1.781387D3*TI)
      RF(13) = 1.D13
      RF(14) = EXP(2.85473118D1 -2.40537567D4*TI)
      RF(15) = EXP(3.22361913D1 -2.01286667D4*TI)
      RF(16) = EXP(4.24761511D1 -8.6D-1*ALOGT)
      TMP = EXP(-1.72D0*ALOGT)
      RF(17) = 3.D20 * TMP
      RF(19) = 3.75D20 * TMP
      RF(18) = EXP(4.36851114D1 -7.6D-1*ALOGT)
      RF(20) = EXP(3.20498617D1 -7.25286183D3*TI)
      RF(21) = 1.D18*TI
      RF(22) = EXP(3.90385861D1 -6.D-1*ALOGT)
      RF(23) = EXP(4.55408762D1 -1.25D0*ALOGT)
      RF(24) = 5.5D20*TI2
      RF(25) = 2.2D22*TI2
      RF(26) = EXP(2.90097872D1 -3.37658384D2*TI)
      RF(27) = EXP(3.09632256D1 -5.37435401D2*TI)
      RF(28) = EXP(3.25288609D1 -3.19542584D2*TI)
      RF(29) = EXP(3.77576522D1 -8.D-1*ALOGT)
      RF(30) = EXP(3.70803784D1 -6.3D-1*ALOGT -1.92731984D2*TI)
      RF(31) = EXP(2.03077504D1 +1.62D0*ALOGT -5.45486868D3*TI)
      RF(32) = EXP(2.77171988D1 +4.8D-1*ALOGT +1.30836334D2*TI)
      RF(33) = 7.34D13
      RF(34) = EXP(2.7014835D1 +4.54D-1*ALOGT -1.81158D3*TI)
      RF(35) = EXP(2.38587601D1 +1.05D0*ALOGT -1.64803459D3*TI)
      RF(36) = 2.D13
      RF(37) = 1.2D13
      RF(38) = 6.D12
      RF(39) = EXP(1.75767107D1 +1.5D0*ALOGT -4.00560467D4*TI)
      RF(40) = EXP(1.9190789D1 +1.51D0*ALOGT -1.72603317D3*TI)
      RF(41) = EXP(1.0482906D1 +2.4D0*ALOGT +1.06178717D3*TI)
      RF(42) = EXP(3.09983169D1 +2.51608334D2*TI)
      RF(43) = 2.D13
      RF(44) = 3.D13
      RF(45) = EXP(1.78408622D1 +1.6D0*ALOGT -2.72743434D3*TI)
      RF(46) = 2.501D13
      RF(47) = EXP(1.84206807D1 +1.6D0*ALOGT -1.570036D3*TI)
      RF(48) = EXP(1.76783433D1 +1.228D0*ALOGT -3.52251667D1*TI)
      RF(49) = 5.D13
      RF(50) = EXP(2.19558261D1 +1.18D0*ALOGT +2.2493785D2*TI)
      RF(51) = 5.D12
      RF(52) = 2.D13
      RF(53) = 1.D12
      RF(54) = EXP(3.26416564D1 -1.18759134D4*TI)
      RF(55) = EXP(3.02112379D1 -7.54825001D2*TI)
      RF(56) = EXP(1.31223634D1 +2.D0*ALOGT -3.63825651D3*TI)
      RF(57) = EXP(1.47156719D1 +2.D0*ALOGT -4.16160184D3*TI)
      RF(58) = EXP(3.03390713D1 -3.01930001D2*TI)
      RF(59) = 2.8D13
      RF(60) = 1.2D13
      RF(61) = 7.D13
      RF(62) = 3.D13
      RF(63) = EXP(3.04036098D1 +2.86833501D2*TI)
      RF(64) = 9.D12
      RF(65) = 7.D12
      RF(66) = 1.4D13
      RF(67) = EXP(2.43067848D1 -4.49875701D3*TI)
      RF(68) = 2.648D13
      RF(69) = EXP(8.10772006D0 +2.81D0*ALOGT -2.94884967D3*TI)
      TMP = EXP(-1.D0*ALOGT -8.55468335D3*TI )
      RF(70) = 2.244D18 * TMP
      RF(71) = 1.87D17 * TMP
      RF(72) = EXP(2.96591694D1 -2.01286667D2*TI)
      RF(73) = EXP(3.05213929D1 -4.52895001D2*TI)
C
C	 thermal data
C
      TN1 = ALOGT - 1.0
      TN2 = T
      TN3 = TN2*T
      TN4 = TN3*T
      TN5 = TN4*T
C
      IF (T .GT. 1.D3) THEN
C
      SMH(1) = -3.20502331D+00 + 9.50158922D+02*TI 
     *         + 3.33727920D+00*TN1 - 2.47012365D-05*TN2 
     *         + 8.32427963D-08*TN3 - 1.49638662D-11*TN4 
     *         + 1.00127688D-15*TN5 
      SMH(2) = -4.46682914D-01 - 2.54736599D+04*TI 
     *         + 2.50000001D+00*TN1 - 1.15421486D-11*TN2 
     *         + 2.69269913D-15*TN3 - 3.94596029D-19*TN4 
     *         + 2.49098679D-23*TN5 
      SMH(3) = 4.78433864D+00 - 2.92175791D+04*TI 
     *         + 2.56942078D+00*TN1 - 4.29870569D-05*TN2 
     *         + 6.99140982D-09*TN3 - 8.34814992D-13*TN4 
     *         + 6.14168455D-17*TN5 
      SMH(4) = 5.45323129D+00 + 1.08845772D+03*TI 
     *         + 3.28253784D+00*TN1 + 7.41543770D-04*TN2 
     *         - 1.26327778D-07*TN3 + 1.74558796D-11*TN4 
     *         - 1.08358897D-15*TN5 
      SMH(5) = 4.47669610D+00 - 3.85865700D+03*TI 
     *         + 3.09288767D+00*TN1 + 2.74214858D-04*TN2 
     *         + 2.10842047D-08*TN3 - 7.32884630D-12*TN4 
     *         + 5.87061880D-16*TN5 
      SMH(6) = 4.96677010D+00 + 3.00042971D+04*TI 
     *         + 3.03399249D+00*TN1 + 1.08845902D-03*TN2 
     *         - 2.73454197D-08*TN3 - 8.08683225D-12*TN4 
     *         + 8.41004960D-16*TN5 
      SMH(7) = 3.78510215D+00 - 1.11856713D+02*TI 
     *         + 4.01721090D+00*TN1 + 1.11991007D-03*TN2 
     *         - 1.05609692D-07*TN3 + 9.52053083D-12*TN4 
     *         - 5.39542675D-16*TN5 
      SMH(8) = 6.17119324D+00 - 4.62636040D+04*TI 
     *         + 2.87410113D+00*TN1 + 1.82819646D-03*TN2 
     *         - 2.34824328D-07*TN3 + 2.16816291D-11*TN4 
     *         - 9.38637835D-16*TN5 
      SMH(9) = 8.62650169D+00 - 5.09259997D+04*TI 
     *         + 2.29203842D+00*TN1 + 2.32794319D-03*TN2 
     *         - 3.35319912D-07*TN3 + 3.48255000D-11*TN4 
     *         - 1.69858183D-15*TN5 
      SMH(10) = 8.48007179D+00 - 1.67755843D+04*TI 
     *         + 2.28571772D+00*TN1 + 3.61995018D-03*TN2 
     *         - 4.97857247D-07*TN3 + 4.96403870D-11*TN4 
     *         - 2.33577197D-15*TN5 
      SMH(11) = 1.84373180D+01 + 9.46834459D+03*TI 
     *         + 7.48514950D-02*TN1 + 6.69547335D-03*TN2 
     *         - 9.55476348D-07*TN3 + 1.01910446D-10*TN4 
     *         - 5.09076150D-15*TN5 
      SMH(12) = 7.81868772D+00 + 1.41518724D+04*TI 
     *         + 2.71518561D+00*TN1 + 1.03126372D-03*TN2 
     *         - 1.66470962D-07*TN3 + 1.91710840D-11*TN4 
     *         - 1.01823858D-15*TN5 
      SMH(13) = 2.27163806D+00 + 4.87591660D+04*TI 
     *         + 3.85746029D+00*TN1 + 2.20718513D-03*TN2 
     *         - 3.69135673D-07*TN3 + 4.36241823D-11*TN4 
     *         - 2.36042082D-15*TN5 
      SMH(14) = 9.79834492D+00 - 4.01191815D+03*TI 
     *         + 2.77217438D+00*TN1 + 2.47847763D-03*TN2 
     *         - 4.14076022D-07*TN3 + 4.90968148D-11*TN4 
     *         - 2.66754356D-15*TN5 
      SMH(15) = 1.36563230D+01 + 1.39958323D+04*TI 
     *         + 1.76069008D+00*TN1 + 4.60000041D-03*TN2 
     *         - 7.37098022D-07*TN3 + 8.38676767D-11*TN4 
     *         - 4.41927820D-15*TN5 
      SMH(16) = 5.81043215D+00 + 3.24250627D+03*TI 
     *         + 3.69266569D+00*TN1 + 4.32288399D-03*TN2 
     *         - 6.25168533D-07*TN3 + 6.56028863D-11*TN4 
     *         - 3.24277101D-15*TN5 
C
      ELSE
C
      SMH(1) = 6.83010238D-01 + 9.17935173D+02*TI 
     *         + 2.34433112D+00*TN1 + 3.99026037D-03*TN2 
     *         - 3.24635850D-06*TN3 + 1.67976745D-09*TN4 
     *         - 3.68805881D-13*TN5 
      SMH(2) = -4.46682853D-01 - 2.54736599D+04*TI 
     *         + 2.50000000D+00*TN1 + 3.52666409D-13*TN2 
     *         - 3.32653273D-16*TN3 + 1.91734693D-19*TN4 
     *         - 4.63866166D-23*TN5 
      SMH(3) = 2.05193346D+00 - 2.91222592D+04*TI 
     *         + 3.16826710D+00*TN1 - 1.63965942D-03*TN2 
     *         + 1.10717733D-06*TN3 - 5.10672187D-10*TN4 
     *         + 1.05632986D-13*TN5 
      SMH(4) = 3.65767573D+00 + 1.06394356D+03*TI 
     *         + 3.78245636D+00*TN1 - 1.49836708D-03*TN2 
     *         + 1.64121700D-06*TN3 - 8.06774591D-10*TN4 
     *         + 1.62186419D-13*TN5 
      SMH(5) = -1.03925458D-01 - 3.61508056D+03*TI 
     *         + 3.99201543D+00*TN1 - 1.20065876D-03*TN2 
     *         + 7.69656402D-07*TN3 - 3.23427778D-10*TN4 
     *         + 6.82057350D-14*TN5 
      SMH(6) = -8.49032208D-01 + 3.02937267D+04*TI 
     *         + 4.19864056D+00*TN1 - 1.01821705D-03*TN2 
     *         + 1.08673369D-06*TN3 - 4.57330885D-10*TN4 
     *         + 8.85989085D-14*TN5 
      SMH(7) = 3.71666245D+00 - 2.94808040D+02*TI 
     *         + 4.30179801D+00*TN1 - 2.37456025D-03*TN2 
     *         + 3.52638152D-06*TN3 - 2.02303245D-09*TN4 
     *         + 4.64612562D-13*TN5 
      SMH(8) = 1.56253185D+00 - 4.60040401D+04*TI 
     *         + 3.76267867D+00*TN1 + 4.84436072D-04*TN2 
     *         + 4.65816402D-07*TN3 - 3.20909294D-10*TN4 
     *         + 8.43708595D-14*TN5 
      SMH(9) = -7.69118967D-01 - 5.04968163D+04*TI 
     *         + 4.19860411D+00*TN1 - 1.18330710D-03*TN2 
     *         + 1.37216037D-06*TN3 - 5.57346651D-10*TN4 
     *         + 9.71573685D-14*TN5 
      SMH(10) = 1.60456433D+00 - 1.64449988D+04*TI 
     *         + 3.67359040D+00*TN1 + 1.00547588D-03*TN2 
     *         + 9.55036427D-07*TN3 - 5.72597854D-10*TN4 
     *         + 1.27192867D-13*TN5 
      SMH(11) = -4.64130376D+00 + 1.02466476D+04*TI 
     *         + 5.14987613D+00*TN1 - 6.83548940D-03*TN2 
     *         + 8.19667665D-06*TN3 - 4.03952522D-09*TN4 
     *         + 8.33469780D-13*TN5 
      SMH(12) = 3.50840928D+00 + 1.43440860D+04*TI 
     *         + 3.57953347D+00*TN1 - 3.05176840D-04*TN2 
     *         + 1.69469055D-07*TN3 + 7.55838237D-11*TN4 
     *         - 4.52212249D-14*TN5 
      SMH(13) = 9.90105222D+00 + 4.83719697D+04*TI 
     *         + 2.35677352D+00*TN1 + 4.49229839D-03*TN2 
     *         - 1.18726045D-06*TN3 + 2.04932518D-10*TN4 
     *         - 7.18497740D-15*TN5 
      SMH(14) = 3.39437243D+00 - 3.83956496D+03*TI 
     *         + 4.22118584D+00*TN1 - 1.62196266D-03*TN2 
     *         + 2.29665743D-06*TN3 - 1.10953411D-09*TN4 
     *         + 2.16884433D-13*TN5 
      SMH(15) = 6.02812900D-01 + 1.43089567D+04*TI 
     *         + 4.79372315D+00*TN1 - 4.95416685D-03*TN2 
     *         + 6.22033347D-06*TN3 - 3.16071051D-09*TN4 
     *         + 6.58863260D-13*TN5 
      SMH(16) = 5.47302243D+00 + 3.19391367D+03*TI 
     *         + 3.86388918D+00*TN1 + 2.79836152D-03*TN2 
     *         + 9.88786318D-07*TN3 - 8.71100100D-10*TN4 
     *         + 2.18483639D-13*TN5 
      ENDIF
C
C	 equilibrium constants 
C
      DO N = 1, 16
          EG(N) = EXP(MIN(SMH(N), 230.26))
      ENDDO
C
      PFAC = PATM/(RU*T)
      PFAC2 = PFAC*PFAC
      PFAC3 = PFAC2*PFAC
C
      EQK(1)=EG(4)/EG(3)/EG(3)/PFAC
      EQK(2)=EG(5)/EG(2)/EG(3)/PFAC
      EQK(3)=EG(2)*EG(5)/EG(1)/EG(3)
      EQK(4)=EG(4)*EG(5)/EG(3)/EG(7)
      EQK(5)=EG(2)*EG(14)/EG(3)/EG(8)
      EQK(6)=EG(1)*EG(12)/EG(3)/EG(9)
      EQK(7)=EG(2)*EG(15)/EG(3)/EG(10)
      EQK(8)=EG(5)*EG(10)/EG(3)/EG(11)
      EQK(9)=EG(13)/EG(3)/EG(12)/PFAC
      EQK(10)=EG(5)*EG(12)/EG(3)/EG(14)
      EQK(11)=EG(2)*EG(13)/EG(3)/EG(14)
      EQK(12)=EG(5)*EG(14)/EG(3)/EG(15)
      EQK(13)=EG(5)*EG(15)/EG(3)/EG(16)
      EQK(14)=EG(3)*EG(13)/EG(4)/EG(12)
      EQK(15)=EG(7)*EG(14)/EG(4)/EG(15)
      EQK(16)=EG(7)/EG(2)/EG(4)/PFAC
      EQK(17)=EQK(16)
      EQK(18)=EQK(16)
      EQK(19)=EQK(16)
      EQK(20)=EG(3)*EG(5)/EG(2)/EG(4)
      EQK(21)=EG(1)/EG(2)/EG(2)/PFAC
      EQK(22)=EQK(21)
      EQK(23)=EQK(21)
      EQK(24)=EQK(21)
      EQK(25)=EG(6)/EG(2)/EG(5)/PFAC
      EQK(26)=EG(3)*EG(6)/EG(2)/EG(7)
      EQK(27)=EG(1)*EG(4)/EG(2)/EG(7)
      EQK(28)=EG(5)*EG(5)/EG(2)/EG(7)
      EQK(29)=EG(10)/EG(2)/EG(8)/PFAC
      EQK(30)=EG(11)/EG(2)/EG(10)/PFAC
      EQK(31)=EG(1)*EG(10)/EG(2)/EG(11)
      EQK(32)=EG(15)/EG(2)/EG(14)/PFAC
      EQK(33)=EG(1)*EG(12)/EG(2)/EG(14)
      EQK(34)=EG(16)/EG(2)/EG(15)/PFAC
      EQK(35)=EG(1)*EG(14)/EG(2)/EG(15)
      EQK(36)=EG(1)*EG(15)/EG(2)/EG(16)
      EQK(37)=EG(5)*EG(10)/EG(2)/EG(16)
      EQK(38)=EG(6)*EG(9)/EG(2)/EG(16)
      EQK(39)=EG(15)/EG(1)/EG(12)/PFAC
      EQK(40)=EG(2)*EG(6)/EG(1)/EG(5)
      EQK(41)=EG(3)*EG(6)/EG(5)/EG(5)
      EQK(42)=EG(4)*EG(6)/EG(5)/EG(7)
      EQK(43)=EG(2)*EG(15)/EG(5)/EG(8)
      EQK(44)=EG(2)*EG(15)/EG(5)/EG(9)
      EQK(45)=EG(6)*EG(8)/EG(5)/EG(10)
      EQK(46)=EG(6)*EG(9)/EG(5)/EG(10)
      EQK(47)=EG(6)*EG(10)/EG(5)/EG(11)
      EQK(48)=EG(2)*EG(13)/EG(5)/EG(12)
      EQK(49)=EG(6)*EG(12)/EG(5)/EG(14)
      EQK(50)=EG(6)*EG(14)/EG(5)/EG(15)
      EQK(51)=EG(6)*EG(15)/EG(5)/EG(16)
      EQK(52)=EG(5)*EG(15)/EG(7)/EG(8)
      EQK(53)=EG(4)*EG(11)/EG(7)/EG(10)
      EQK(54)=EG(5)*EG(13)/EG(7)/EG(12)
      EQK(55)=EG(5)*EG(14)/EG(4)/EG(8)
      EQK(56)=EG(2)*EG(10)/EG(1)/EG(8)
      EQK(57)=EG(10)*EG(10)/EG(8)/EG(11)
      EQK(58)=EG(8)/EG(9)
      EQK(62)=EQK(58)
      EQK(64)=EQK(58)
      EQK(65)=EQK(58)
      EQK(59)=EG(2)*EG(5)*EG(12)/EG(4)/EG(9)*PFAC
      EQK(60)=EG(6)*EG(12)/EG(4)/EG(9)
      EQK(61)=EG(2)*EG(10)/EG(1)/EG(9)
      EQK(63)=EG(10)*EG(10)/EG(9)/EG(11)
      EQK(66)=EG(12)*EG(15)/EG(9)/EG(13)
      EQK(67)=EG(5)*EG(15)/EG(4)/EG(10)
      EQK(68)=EG(11)*EG(12)/EG(10)/EG(14)
      EQK(69)=EG(11)*EG(14)/EG(10)/EG(15)
      EQK(70)=EG(2)*EG(12)/EG(14)*PFAC
      EQK(71)=EQK(70)
      EQK(72)=EG(7)*EG(12)/EG(4)/EG(14)
      EQK(73)=EG(7)*EG(15)/EG(4)/EG(16)
C
C	Compute reverse reaction rates
C
      DO I = 1, 73
        RB(I) = RF(I) / MAX(EQK(I),1.0D-100)
      ENDDO
C
C	 rates at low pressure limit
C
      RKLOW(1) = EXP(6.33329483D1 -3.14D0*ALOGT -6.18956501D2/T)
      RKLOW(2) = EXP(7.68923562D1 -4.76D0*ALOGT -1.22784867D3/T)
      RKLOW(3) = EXP(5.55621468D1 -2.57D0*ALOGT -7.17083751D2/T)
      RKLOW(4) = EXP(7.39217399D1 -4.82D0*ALOGT -3.28600484D3/T)
      RKLOW(5) = EXP(6.37931383D1 -3.42D0*ALOGT -4.24463259D4/T)
C
C      third-body concentrations
C
      CTOT = 0.0
      DO K = 1, 13
         CTOT = CTOT + C(K)
      ENDDO
C
      CTB(1) = CTOT + 1.4D0*C(1) + 1.44D1*C(6) 
     *         + C(9) + 7.5D-1*C(10) + 2.6D0*C(11) 
      CTB(2) = CTOT + C(1) + 5.D0*C(6) 
     *         + C(9) + 5.D-1*C(10) + C(11) 
      CTB(9) = CTOT + C(1) + 5.D0*C(4) 
     *         + 5.D0*C(6) + C(9) + 5.D-1*C(10) 
     *         + 2.5D0*C(11) 
      CTB(16) = CTOT - C(4) - C(6) 
     *         - 2.5D-1*C(10) + 5.D-1*C(11) - C(13) 
      CTB(21) = CTOT - C(1) - C(6) 
     *         + C(9) - C(11) 
      CTB(25) = CTOT - 2.7D-1*C(1) + 2.65D0*C(6) 
     *         + C(9) 
      CTB(29) = CTOT + C(1) + 5.D0*C(6) 
     *         + C(9) + 5.D-1*C(10) + C(11) 
      CTB(30) = CTB(29) 
      CTB(32) = CTB(29) 
      CTB(34) = CTB(29) 
      CTB(39) = CTB(29) 
      CTB(71) = CTOT + C(1) - C(6) 
     *         + C(9) + 5.D-1*C(10) + C(11) 
C
C      fall-off reactions
C
C
      PR = RKLOW(1) * CTB(29) / RF(29)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 3.2D-1*EXP(-T/7.8D1) + 6.8D-1*EXP(-T/1.995D3)
     *     + EXP(-5.59D3/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(29) = RF(29) * PCOR
      RB(29) = RB(29) * PCOR
C
      PR = RKLOW(2) * CTB(30) / RF(30)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 2.17D-1*EXP(-T/7.4D1) + 7.83D-1*EXP(-T/2.941D3)
     *     + EXP(-6.964D3/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(30) = RF(30) * PCOR
      RB(30) = RB(30) * PCOR
C
      PR = RKLOW(3) * CTB(32) / RF(32)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 2.176D-1*EXP(-T/2.71D2) + 7.824D-1*EXP(-T/2.755D3)
     *     + EXP(-6.57D3/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(32) = RF(32) * PCOR
      RB(32) = RB(32) * PCOR
C
      PR = RKLOW(4) * CTB(34) / RF(34)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 2.813D-1*EXP(-T/1.03D2) + 7.187D-1*EXP(-T/1.291D3)
     *     + EXP(-4.16D3/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(34) = RF(34) * PCOR
      RB(34) = RB(34) * PCOR
C
      PR = RKLOW(5) * CTB(39) / RF(39)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 6.8D-2*EXP(-T/1.97D2) + 9.32D-1*EXP(-T/1.54D3)
     *     + EXP(-1.03D4/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(39) = RF(39) * PCOR
      RB(39) = RB(39) * PCOR
C
      RF(1) = RF(1)*CTB(1)*C(3)*C(3)
      RF(2) = RF(2)*CTB(2)*C(3)*C(2)
      RF(3) = RF(3)*C(3)*C(1)
      RF(4) = RF(4)*C(3)*C(7)
      RF(5) = RF(5)*C(3)
      RF(6) = RF(6)*C(3)
      RF(7) = RF(7)*C(3)*C(8)
      RF(8) = RF(8)*C(3)*C(9)
      RF(9) = RF(9)*CTB(9)*C(3)*C(10)
      RF(10) = RF(10)*C(3)
      RF(11) = RF(11)*C(3)
      RF(12) = RF(12)*C(3)*C(12)
      RF(13) = RF(13)*C(3)
      RF(14) = RF(14)*C(4)*C(10)
      RF(15) = RF(15)*C(4)*C(12)
      RF(16) = RF(16)*CTB(16)*C(2)*C(4)
      RF(17) = RF(17)*C(2)*C(4)*C(4)
      RF(18) = RF(18)*C(2)*C(4)*C(6)
      RF(19) = RF(19)*C(2)*C(4)*C(13)
      RF(20) = RF(20)*C(2)*C(4)
      RF(21) = RF(21)*CTB(21)*C(2)*C(2)
      RF(22) = RF(22)*C(2)*C(2)*C(1)
      RF(23) = RF(23)*C(2)*C(2)*C(6)
      RF(24) = RF(24)*C(2)*C(2)*C(11)
      RF(25) = RF(25)*CTB(25)*C(2)*C(5)
      RF(26) = RF(26)*C(2)*C(7)
      RF(27) = RF(27)*C(2)*C(7)
      RF(28) = RF(28)*C(2)*C(7)
      RF(29) = RF(29)*C(2)
      RF(30) = RF(30)*C(2)*C(8)
      RF(31) = RF(31)*C(2)*C(9)
      RF(32) = RF(32)*C(2)
      RF(33) = RF(33)*C(2)
      RF(34) = RF(34)*C(2)*C(12)
      RF(35) = RF(35)*C(2)*C(12)
      RF(36) = RF(36)*C(2)
      RF(37) = RF(37)*C(2)
      RF(38) = RF(38)*C(2)
      RF(39) = RF(39)*C(1)*C(10)
      RF(40) = RF(40)*C(5)*C(1)
      RF(41) = RF(41)*C(5)*C(5)
      RF(42) = RF(42)*C(5)*C(7)
      RF(43) = RF(43)*C(5)
      RF(44) = RF(44)*C(5)
      RF(45) = RF(45)*C(5)*C(8)
      RF(46) = RF(46)*C(5)*C(8)
      RF(47) = RF(47)*C(5)*C(9)
      RF(48) = RF(48)*C(5)*C(10)
      RF(49) = RF(49)*C(5)
      RF(50) = RF(50)*C(5)*C(12)
      RF(51) = RF(51)*C(5)
      RF(52) = RF(52)*C(7)
      RF(53) = RF(53)*C(7)*C(8)
      RF(54) = RF(54)*C(7)*C(10)
      RF(55) = RF(55)*C(4)
      RF(56) = RF(56)*C(1)
      RF(57) = RF(57)*C(9)
      RF(58) = RF(58)*C(13)
      RF(59) = RF(59)*C(4)
      RF(60) = RF(60)*C(4)
      RF(61) = RF(61)*C(1)
      RF(62) = RF(62)*C(6)
      RF(63) = RF(63)*C(9)
      RF(64) = RF(64)*C(10)
      RF(65) = RF(65)*C(11)
      RF(66) = RF(66)*C(11)
      RF(67) = RF(67)*C(8)*C(4)
      RF(68) = RF(68)*C(8)
      RF(69) = RF(69)*C(8)*C(12)
      RF(70) = RF(70)*C(6)
      RF(71) = RF(71)*CTB(71)
      RF(72) = RF(72)*C(4)
      RF(73) = RF(73)*C(4)
      RB(1) = RB(1)*CTB(1)*C(4)
      RB(2) = RB(2)*CTB(2)*C(5)
      RB(3) = RB(3)*C(2)*C(5)
      RB(4) = RB(4)*C(5)*C(4)
      RB(5) = RB(5)*C(2)
      RB(6) = RB(6)*C(1)*C(10)
      RB(7) = RB(7)*C(2)*C(12)
      RB(8) = RB(8)*C(5)*C(8)
      RB(9) = RB(9)*CTB(9)*C(11)
      RB(10) = RB(10)*C(5)*C(10)
      RB(11) = RB(11)*C(2)*C(11)
      RB(12) = RB(12)*C(5)
      RB(13) = RB(13)*C(5)*C(12)
      RB(14) = RB(14)*C(3)*C(11)
      RB(15) = RB(15)*C(7)
      RB(16) = RB(16)*CTB(16)*C(7)
      RB(17) = RB(17)*C(7)*C(4)
      RB(18) = RB(18)*C(7)*C(6)
      RB(19) = RB(19)*C(7)*C(13)
      RB(20) = RB(20)*C(3)*C(5)
      RB(21) = RB(21)*CTB(21)*C(1)
      RB(22) = RB(22)*C(1)*C(1)
      RB(23) = RB(23)*C(1)*C(6)
      RB(24) = RB(24)*C(1)*C(11)
      RB(25) = RB(25)*CTB(25)*C(6)
      RB(26) = RB(26)*C(3)*C(6)
      RB(27) = RB(27)*C(4)*C(1)
      RB(28) = RB(28)*C(5)*C(5)
      RB(29) = RB(29)*C(8)
      RB(30) = RB(30)*C(9)
      RB(31) = RB(31)*C(8)*C(1)
      RB(32) = RB(32)*C(12)
      RB(33) = RB(33)*C(1)*C(10)
      RB(35) = RB(35)*C(1)
      RB(36) = RB(36)*C(1)*C(12)
      RB(37) = RB(37)*C(5)*C(8)
      RB(38) = RB(38)*C(6)
      RB(39) = RB(39)*C(12)
      RB(40) = RB(40)*C(2)*C(6)
      RB(41) = RB(41)*C(3)*C(6)
      RB(42) = RB(42)*C(4)*C(6)
      RB(43) = RB(43)*C(2)*C(12)
      RB(44) = RB(44)*C(2)*C(12)
      RB(45) = RB(45)*C(6)
      RB(46) = RB(46)*C(6)
      RB(47) = RB(47)*C(8)*C(6)
      RB(48) = RB(48)*C(2)*C(11)
      RB(49) = RB(49)*C(6)*C(10)
      RB(50) = RB(50)*C(6)
      RB(51) = RB(51)*C(6)*C(12)
      RB(52) = RB(52)*C(5)*C(12)
      RB(53) = RB(53)*C(4)*C(9)
      RB(54) = RB(54)*C(5)*C(11)
      RB(55) = RB(55)*C(5)
      RB(56) = RB(56)*C(2)*C(8)
      RB(57) = RB(57)*C(8)*C(8)
      RB(58) = RB(58)*C(13)
      RB(59) = RB(59)*C(2)*C(5)*C(10)
      RB(60) = RB(60)*C(10)*C(6)
      RB(61) = RB(61)*C(8)*C(2)
      RB(62) = RB(62)*C(6)
      RB(63) = RB(63)*C(8)*C(8)
      RB(64) = RB(64)*C(10)
      RB(65) = RB(65)*C(11)
      RB(66) = RB(66)*C(10)*C(12)
      RB(67) = RB(67)*C(5)*C(12)
      RB(68) = RB(68)*C(9)*C(10)
      RB(69) = RB(69)*C(9)
      RB(70) = RB(70)*C(2)*C(10)*C(6)
      RB(71) = RB(71)*CTB(71)*C(2)*C(10)
      RB(72) = RB(72)*C(7)*C(10)
      RB(73) = RB(73)*C(7)*C(12)
C
C	 solving QSS species concentration
C	-------------------------------------------
C     CH2
      ABV1 = RB(29) + RB(43) + RF(45)
     *          + RB(52) + RB(56) + RB(57)
      DEN1 =  RF(5) +  RF(29) +  RF(43)
     *          +  RB(45) +  RF(52) +  RF(55)
     *          +  RF(56) +  RF(57) +  RB(58)
     *          +  RB(62) +  RB(64) +  RB(65)
C
C     CH2(S)
      ABV2 = RB(6) + RB(44) + RF(46)
     *          + RB(59) + RB(60) + RB(61)
     *          + RB(63) + RB(66)
      DEN2 =  RF(6) +  RB(38) +  RF(44)
     *          +  RB(46) +  RF(58) +  RF(59)
     *          +  RF(60) +  RF(61) +  RF(62)
     *          +  RF(63) +  RF(64) +  RF(65)
     *          +  RF(66)
C
C     HCO
      ABV3 = RB(10) + RB(11) + RF(12)
     *          + RF(15) + RB(32) + RB(33)
     *          + RF(35) + RB(49) + RF(50)
     *          + RB(68) + RF(69) + RB(70)
     *          + RB(71) + RB(72)
      DEN3 =  RB(5) +  RF(10) +  RF(11)
     *          +  RB(12) +  RB(15) +  RF(32)
     *          +  RF(33) +  RB(35) +  RF(49)
     *          +  RB(50) +  RB(55) +  RF(68)
     *          +  RB(69) +  RF(70) +  RF(71)
     *          +  RF(72)
C
C     CH2OH
      ABV4 = RB(13) + RF(34) + RB(36)
     *          + RB(37) + RB(51) + RB(73)
      DEN4 =  RF(13) +  RB(34) +  RF(36)
     *          +  RF(37) +  RF(38) +  RF(51)
     *          +  RF(73)
C
	F1 = ABV1/DEN1
	A12 = (RF(58)+RF(62)+RF(64)+RF(65))/DEN1
	A13 = (RB(5)+RB(55))/DEN1
C
	F2 = ABV2/DEN2
	A21 = (RB(58)+RB(62)+RB(64)+RB(65))/DEN2
	A24 = RF(38)/DEN2
C
	F3 = ABV3/DEN3
	A31 = (RF(5)+RF(55))/DEN3
C
	F4 = ABV4/DEN4
	A42 = RB(38)/DEN4
C
	TEMP = 1 - A13*A31
	FF1 = (F1 + A13*F3)/TEMP
	AA12 = A12/TEMP
C
	TEMP = 1 - A24*A42
	FF2 = (F2 + A24*F4)/TEMP
	AA21 = A21/TEMP
C
	XQ(1) = (FF1 + AA12*FF2)/(1 - AA12*AA21)
	XQ(2) = FF2 + AA21*XQ(1)
	XQ(3) = F3 + A31*XQ(1)
	XQ(4) = F4 + A42*XQ(2)
C
C	---------------------------------------------
C
C	 update rates of reactions involving QSS species
C
      RF(5) = RF(5) * XQ(1)
      RB(5) = RB(5) * XQ(3)
      RF(6) = RF(6) * XQ(2)
      RF(10) = RF(10) * XQ(3)
      RF(11) = RF(11) * XQ(3)
      RB(12) = RB(12) * XQ(3)
      RF(13) = RF(13) * XQ(4)
      RB(15) = RB(15) * XQ(3)
      RF(29) = RF(29) * XQ(1)
      RF(32) = RF(32) * XQ(3)
      RF(33) = RF(33) * XQ(3)
      RB(34) = RB(34) * XQ(4)
      RB(35) = RB(35) * XQ(3)
      RF(36) = RF(36) * XQ(4)
      RF(37) = RF(37) * XQ(4)
      RF(38) = RF(38) * XQ(4)
      RB(38) = RB(38) * XQ(2)
      RF(43) = RF(43) * XQ(1)
      RF(44) = RF(44) * XQ(2)
      RB(45) = RB(45) * XQ(1)
      RB(46) = RB(46) * XQ(2)
      RF(49) = RF(49) * XQ(3)
      RB(50) = RB(50) * XQ(3)
      RF(51) = RF(51) * XQ(4)
      RF(52) = RF(52) * XQ(1)
      RF(55) = RF(55) * XQ(1)
      RB(55) = RB(55) * XQ(3)
      RF(56) = RF(56) * XQ(1)
      RF(57) = RF(57) * XQ(1)
      RF(58) = RF(58) * XQ(2)
      RB(58) = RB(58) * XQ(1)
      RF(59) = RF(59) * XQ(2)
      RF(60) = RF(60) * XQ(2)
      RF(61) = RF(61) * XQ(2)
      RF(62) = RF(62) * XQ(2)
      RB(62) = RB(62) * XQ(1)
      RF(63) = RF(63) * XQ(2)
      RF(64) = RF(64) * XQ(2)
      RB(64) = RB(64) * XQ(1)
      RF(65) = RF(65) * XQ(2)
      RB(65) = RB(65) * XQ(1)
      RF(66) = RF(66) * XQ(2)
      RF(68) = RF(68) * XQ(3)
      RB(69) = RB(69) * XQ(3)
      RF(70) = RF(70) * XQ(3)
      RF(71) = RF(71) * XQ(3)
      RF(72) = RF(72) * XQ(3)
      RF(73) = RF(73) * XQ(4)
C
C	Compute WDOT
C
      DO I = 1, 73
         ROP(I) = RF(I) - RB(I)
      ENDDO
C
C	Combine duplicated reactions
C
      ROP(16) = ROP(16) + ROP(17) + ROP(18) + ROP(19)
	ROP(21) = ROP(21) + ROP(22) + ROP(23) + ROP(24)
      ROP(58) = ROP(58) + ROP(62) + ROP(64) + ROP(65)
      ROP(70) = ROP(70) + ROP(71)
C
      WDOT(1) = -ROP(3) +ROP(6) +ROP(21) +ROP(27) 
     *          +ROP(31) +ROP(33) +ROP(35) +ROP(36) 
     *          -ROP(39) -ROP(40) -ROP(56) -ROP(61) 
      WDOT(2) = -ROP(2) +ROP(3) +ROP(5) +ROP(7) 
     *          +ROP(11) -ROP(16) -ROP(20) -ROP(21) -ROP(21) 
     *          -ROP(25) -ROP(26) -ROP(27) -ROP(28) 
     *          -ROP(29) -ROP(30) -ROP(31) -ROP(32) 
     *          -ROP(33) -ROP(34) -ROP(35) -ROP(36) 
     *          -ROP(37) -ROP(38) +ROP(40) +ROP(43) 
     *          +ROP(44) +ROP(48) +ROP(56) +ROP(59) 
     *          +ROP(61) +ROP(70) 
      WDOT(3) = -ROP(1) -ROP(1) -ROP(2) -ROP(3) 
     *          -ROP(4) -ROP(5) -ROP(6) -ROP(7) 
     *          -ROP(8) -ROP(9) -ROP(10) -ROP(11) 
     *          -ROP(12) -ROP(13) +ROP(14) +ROP(20) 
     *          +ROP(26) +ROP(41) 
      WDOT(4) = +ROP(1) +ROP(4) -ROP(14) -ROP(15) 
     *          -ROP(16) -ROP(20) +ROP(27) +ROP(42) 
     *          +ROP(53) -ROP(55) -ROP(59) -ROP(60) 
     *          -ROP(67) -ROP(72) -ROP(73) 
      WDOT(5) = +ROP(2) +ROP(3) +ROP(4) +ROP(8) 
     *          +ROP(10) +ROP(12) +ROP(13) +ROP(20) 
     *          -ROP(25) +ROP(28) +ROP(28) +ROP(37) 
     *          -ROP(40) -ROP(41) -ROP(41) -ROP(42) 
     *          -ROP(43) -ROP(44) -ROP(45) -ROP(46) 
     *          -ROP(47) -ROP(48) -ROP(49) -ROP(50) 
     *          -ROP(51) +ROP(52) +ROP(54) +ROP(55) 
     *          +ROP(59) +ROP(67) 
      WDOT(6) = +ROP(25) +ROP(26) +ROP(38) +ROP(40) 
     *          +ROP(41) +ROP(42) +ROP(45) +ROP(46) 
     *          +ROP(47) +ROP(49) +ROP(50) +ROP(51) 
     *          +ROP(60) 
      WDOT(7) = -ROP(4) +ROP(15) +ROP(16) -ROP(26) 
     *          -ROP(27) -ROP(28) -ROP(42) -ROP(52) 
     *          -ROP(53) -ROP(54) +ROP(72) +ROP(73) 
      WDOT(8) = -ROP(7) +ROP(8) +ROP(29) -ROP(30) 
     *          +ROP(31) +ROP(37) -ROP(45) -ROP(46) 
     *          +ROP(47) -ROP(53) +ROP(56) +ROP(57) +ROP(57) 
     *          +ROP(61) +ROP(63) +ROP(63) -ROP(67) 
     *          -ROP(68) -ROP(69) 
      WDOT(9) = -ROP(8) +ROP(30) -ROP(31) -ROP(47) 
     *          +ROP(53) -ROP(57) -ROP(63) +ROP(68) 
     *          +ROP(69) 
      WDOT(10) = +ROP(6) -ROP(9) +ROP(10) -ROP(14) 
     *          +ROP(33) -ROP(39) -ROP(48) +ROP(49) 
     *          -ROP(54) +ROP(59) +ROP(60) +ROP(66) 
     *          +ROP(68) +ROP(70) +ROP(72) 
      WDOT(11) = +ROP(9) +ROP(11) +ROP(14) +ROP(48) 
     *          +ROP(54) -ROP(66) 
      WDOT(12) = +ROP(7) -ROP(12) +ROP(13) -ROP(15) 
     *          +ROP(32) -ROP(34) -ROP(35) +ROP(36) 
     *          +ROP(39) +ROP(43) +ROP(44) -ROP(50) 
     *          +ROP(51) +ROP(52) +ROP(66) +ROP(67) 
     *          -ROP(69) +ROP(73) 
      WDOT(13) = 0.0
C
      RETURN
      END
