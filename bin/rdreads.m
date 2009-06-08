 # comp zones only
 # if((k>=0)&&(k<$nz-4)&&(j>=0)&&(j<$ny-4)&&(i>=0)&&(i<$nx-4))
 # all zones
 # if((k>=-2)&&(k<$nz-2)&&(j>=-2)&&(j<$ny-2)&&(i>=-2)&&(i<$nx-2))
 # cut through E0,2,8,10
 # E8,E10
 # if((k>=5)&&(k<$nz-2)&&(j>=-2)&&(j<$ny-2)&&(i>=5)&&(i<$nx-7))
 # E0,E2
 # if((k>=-2)&&(k<$nz-7)&&(j>=-2)&&(j<$ny-2)&&(i>=5)&&(i<$nx-7))
 # cut through E4,5,6,7
 # if((k>=5)&&(k<$nz-7)&&(j>=-2)&&(j<$ny-2)&&(i>=-2)&&(i<$nx-2))
 #
 # if((k>=-2)&&(k<$nz-2)&&(j==$ny/3)&&(i>=-2)&&(i<$nx-2))
 # cut through E3,1,11,9
 # if((k>=-2)&&(k<$nz-2)&&(j>=7)&&(j<=$ny-9)&&(i>=-2)&&(i<$nx-2))
 #
 # if((k==-2)&&(j==$ny/3)&&(i>=-2)&&(i<$nx-2))
 # see other macros: ~/sm/smstart
 #                   /usr/local/lib/sm/macro
 # typical start call:
 # re twod.m gpa 0 rdraft
 # if using older no-ver/type files, then specify type/version manually.
 # define grandpa 1
 # e.g.:
 # define fileversion 1
 # define filetype 3
 #
 #// file versions numbers(use sm for backwards compat)
 #define PVER 6
 #define GRIDVER 2
 #define DVER 1    // dumps same as for pdumps, adumps
 #define FLVER 2
 #define NPVER 2
 #define AVG1DVER 2
 #define AVG2DVER 2
 #define ENERVER 5
 #define LOSSVER 5
 #define SPVER   1
 #define TSVER   1
 #define LOGDTVER 1
 #define STEPVER 1
 #define PERFVER 3
 #define ADVER DVER
 #define PDVER DVER
 #define CALCVER 1
 #// type designations for sm automagical read in correct format for similar things
 #define PTYPE     1 // global par file
 #define GRIDTYPE  2
 #define DTYPE     3 // dump
 #define FLTYPE    4 // floor
 #define NPTYPE    5 // np
 #define AVG2DTYPE 6
 #define AVG1DTYPE 7
 #define ENERTYPE  8
 #define LOSSTYPE  9
 #define SPTYPE    10
 #define TSTYPE    11
 #define LOGDTTYPE 12 
 #define STEPTYPE  13
 #define PERFTYPE  14
 #define ADTYPE    15 // same as dump except filename
 #define PDTYPE    16 // same as dump except filename
 #define CALCTYPE  17 // arbitrary calcs during pp
 #define EXPANDTYPE 50 // used to signify doing pp expansion
 #define NPCOMPUTETYPE 33 // used to signify want to compute np before output
 #
 # ctype default rd dump0001.dat 0 rd adump0001.dat 2 set vxdiff=v0x-v2x pllim 0 x1 vxdiff 150 1100 -1E-10 1E-10
 # ctype red rd dump0001.dat 0 rd adump0001.dat 2 set vxdiff=v0x-v2x plo 0 x1 vxdiff
 # ctype blue rd dump0001.dat 0 rd adump0001.dat 2 set vxdiff=v0x-v2x plo 0 x1 vxdiff
 # ctype green rd dump0001.dat 0 rd adump0001.dat 2 set vxdiff=v0x-v2x plo 0 x1 vxdiff
 # ctype cyan rd dump0001.dat 0 rd adump0001.dat 2 set vxdiff=v0x-v2x plo 0 x1 vxdiff
                   #xla $2
                   #yla $3
                   # relocate (15000 31500)
                   # cp *.m /us1/jon/ndata25
                   # rd /us1/jon/ndata25/0_loss.dat
                   # der t min1 tdn25 min1dn25
                   # ctype blue plo null tdn25 min1dn25
                   # ctype red plo null tdn16 (mx3in1dn16/min1dn16)
                   # relocate (15000 31500)
                   # ctype default label n17 n16 n25 \dot{M}(t)
                   # device postencap /us1/jon/compare/dot/n17n16n25mdot.eps
                   # device postscript
                   # define interp (0) set _gam=(5/3) define gam (_gam) set wgam=1
                   # label n17 n16 n25 \dot{M}(t)
                   #label Mass Accretion rate vs. time
                   #xla tc^3/GM
                   #xla freq GM/c^3
  		   #yla \dot{M}/\dot{M}_{inj}
                   #yla (M/M_{inj})^2
		   #yla \dot{M}_{in} \rho_0 R_0^3 c^3/GM
                   # cp *.m /us1/jon/f8
                   # cat /us1/jon/f8/0_numdumps.dat
                   # device postencap /us1/jon/compare/n28rho.eps
                   # device postscript
                   #
                   #label Power Spectrum of \dot{M_{in}}/\dot{M_{inj}}
                   #label Power Spectrum of \dot{M_{in}}
                   #xla t*GM/c^3
                   #xla Frequency*c^3/GM
                   #yla \dot{M_{in}}/\dot{M_{inj}}
                   #yla Power*(c^2/GM)^2
                   #yla \dot{M_{in}} \rho_0 R_0^3
                   #xla $2
                   #yla $3
                   #
                #lweight 2
                #limits 0 23.1 -23.1 23.1
                #xla Rc^2/GM
                #yla zc^2/GM
                #yla \theta
                #
		#read {x 1 y 2 r 3}
                # ticksize -1 0 0 0
                # device postencap /us1/jon/compare/a1all/1d/r1d1.eps
                # pl 0 x2 en1d1
                # pl 0 (LG(x1)) en1d1
                # relocate (10383 31500) label \rho(\theta)
                # device X11
                #
                # global test functions for injection runs:
                # derived from r=1.05rg..20rg injection runs on alphadog/wiseguy
                # set testrho=exp(-(x2-3.14159*0.5)**2/(0.4**2))*0.012*(x1-2)**(.4)
                # K=en/rho^gam
                # set test1ok=exp(-(x2-3.14159*0.5)**2/(0.5**2))*2.33*(x1-2)**(.51)
                # solve this for en
                # set testen=(1/test1ok)*testrho**($gam)
                # set testvx=-(-.56*LG(x1-2)+.87)**(1/0.25)

                # set testvy=.6*sin(((x2-0.2)/2.75)*2*3.14159)*(-0.4*(LG((x1-2)/(17))))
                # set testvz=(-.6*LG(x1-2)+1.06)*exp(-(x2-3.14159*0.5)**2/(0.7**2))
                
                # set testvz = 1.2512*(x1-2)**(-.571305)*exp(-(x2-3.14159*0.5)**2/(0.7**2))


                #device postencap god.eps
                #device ppm filename.ppm
gpa       1     #
                define grandpa $1
                if($grandpa){ echo grandpa file }
                #
labeltime 0     #
                if(!($finaldraft)){\
                 if( ($filetype==$DTYPE)||($filetype==$ADTYPE)||($filetype==$PDTYPE)||($filetype==$FLTYPE)||($filetype==$NPTYPE)){\
                   relocate (10383 31500)
                   set tempt=sprintf('%5.2g',$time)
                   define temptime (tempt)
                   #label t=$time
                   label t=$temptime
	 	 }
		}
labelaxes 1     # (e.g. labelaxes 0 will just label and not try to define)
		if($1){ defineaxes }
                xla $x1label
                yla $x2label
defineaxes 0    #
		if(('$paxes1'=='x1')||('$paxes1'=='x12')){\
		 if($coord==3){\
		  if($interp==0){\
		   define x1label "r c^2/GM"
		  }
		  if($interp==1){\
		   define x1label "R c^2/GM"
		  }
		 }
		 if($coord==2){\
		  if($interp==0){\
		   define x1label "r"
		  }
		  if($interp==1){\
		   define x1label "r"
		  }
		 }
		 if($coord==1){\
		  if($interp==0){\
		   define x1label "x"
		  }
		  if($interp==1){\
		   define x1label "x"
		  }
		 }
		}
		if(('$paxes2'=='x1')||('$paxes2'=='x12')){\
		 if($interp==0){\
		  define x2label "r c^2/GM"
		 }
		 if($interp==1){\
		  define x2label "R c^2/GM"
		 }
		}
		if(('$paxes1'=='x2')||('$paxes1'=='x22')){\
		 if($interp==0){\
		  define x1label "\theta"
		 }
		 if($interp==1){\
		  define x1label "z c^2/GM"
		 }
		}
		if(('$paxes2'=='x2')||('$paxes2'=='x22')){\
		 if($coord==3){\
		  if($interp==0){\
		   define x2label "\theta"
		  }
		  if($interp==1){\
		   define x2label "z c^2/GM"
		  }
		 }
		 if($coord==2){\
		  if($interp==0){\
		   define x2label "z"
		  }
		  if($interp==1){\
		   define x2label "z"
		  }
		 }
		 if($coord==1){\
		  if($interp==0){\
		   define x2label "y"
		  }
		  if($interp==1){\
		   define x2label "y"
		  }
		 }
		}
                if('$paxes1'=='r'){\
                    define x1label "\rho"
                }
                if('$paxes1'=='freq'){\
                    define x1label "Hz GM/c^3"
                }
                if('$paxes2'=='min1dfpow'){\
		       define x2label "Power (\dot{M}_{inj} GM/c^3)^2"
                }
                if('$paxes2'=='r'){\
                    define x2label "\rho"
                }
                if('$paxes3'=='r'){\
                    define x3label "\rho"
                }
                if( ($filetype==$MODETYPE)||($filetype==$ENERTYPE)||($filetype==$LOSSTYPE)||($filetype==$SPTYPE)||($filetype==$TSTYPE)){\
		 if( ('$paxes1'=='t')||('$paxes1'=='td') ){\
		  define x1label "t c^3/GM"
		 }
		 if('$paxes2'=='min1d'){\
		        #define x2label "\dot{M}/\dot{M}_{inj}"
		        define x2label "\dot{M} c^3/(\rho_0 (GM)^2)"
		 }
		 if('$paxes2'=='(ein1d/min1d)'){\
		        define x2label "\dot{E}^{tot}/\dot{M}"
		 }
		 if('$paxes2'=='(mx3in1d/min1d)'){\
		  define x2label "\dot{L}/\dot{M}"
		 }
		}                
                #
labelaxes3 3    #
                defineaxes
                label3 x $1 $x1label
                label3 y $2 $x2label
                label3 z $3 $x3label
                #
prepaxes 13     # Used since processing can use odd names sent to plot macro
                define paxes1 $1
                define paxes2 $2
                if($?3 == 1){\
                 define paxes3 $3
                }\
                else{\
                 define paxes3 (0)
                }
                #
mybox2d 0       # assumes xl,xh,yl,yh already set(hack for plc since contour doesn't know about grid!)
                # not really correct for anything but pure log grids(lg(r)), i.e. can't do axes right for lg(r-2)
                if( ($interp==1)||($coord==1)||($coord==2) ){\
                 ticksize 0 0 0 0
		 box
	        }\
                else{\
		 if(1==1){\
		  ticksize -1 0 0 0
                  # redo limits for log10 in x1
		  set _newxl=LG($xl)
	      	  define newxl (_newxl)
       		  set _newxh=LG($xh)
	          define newxh (_newxh)
		  set _newyl=$yl
		  define newyl (_newyl)
		  set _newyh=$yh
		  define newyh (_newyh)
		  limits  $newxl $newxh $newyl $newyh
		  box
		 }\
		 else{\
	          ticksize -1 0 -1 0
                  # redo limits for log10 in x1,x2
		  set _newxl=LG($xl)
		  define newxl (_newxl)
		  set _newxh=LG($xh)
		  define newxh (_newxh)
		  set _newyl=LG($yl)
		  define newyl (_newyl)
		  set _newyh=LG($yh)
		  define newyh (_newyh)
		  limits  $newxl $newxh $newyl $newyh
		  box
		 }
                 relocate (15000 31500)
                 label RADIAL SCALE NOT EXACTLY CORRECT
	        }
                #
mybox1d 0       #
                #
                box
                #
mylines    0    #
                define LSTART ($LFINAL+1)
                define LFINAL ($LSTART+1)
                lines $LSTART $LFINAL
                #
filetest    1   # use fact that all non-typed files have more than 2 columns for 2nd line
                # only works back to but not including o11 for dumps, rest is ok
                da $1
                define LSTART (1)
                define LFINAL (2)
                lines $LSTART $LFINAL
                read !{test1 1 test2 2 test3 3}
		if(test3<1.001e+36){\
		 gpa 1
                 # put here which versions likely for old non-versioned files
		 if( ('$1'=='0_gparam.par')||('$1'=='0_igparam.par') ){\
                  # old version is likely 4
 		  define fileversion 4
		  define filetype $PTYPE
                 }
		 if( ('$1'=='0_grid1.par')||('$1'=='0_grid2.par')||('$1'=='0_gridact1.par')||('$1'=='0_gridact2.par')||('$1'=='0_igridact1.par')||('$1'=='0_igridact2.par') ){\
 		  define fileversion 2
		  define filetype $GRIDTYPE
                 }
                 set temptemp=substr('$1',0,4)
		 if((temptemp=='dump')||(temptemp=='idum')){\
 		  define fileversion 1
		  define filetype $DTYPE
                 }
		 if((temptemp=='adum')||(temptemp=='iadu')){\
 		  define fileversion 1
		  define filetype $ADTYPE
                 }
		 if((temptemp=='pdum')||(temptemp=='ipdu')){\
 		  define fileversion 1
		  define filetype $PDTYPE
                 }
		 if((temptemp=='floo')||(temptemp=='iflo')){\
 		  define fileversion 2
		  define filetype $FLTYPE
                 }
		 if((temptemp=='npdu')||(temptemp=='inpd')){\
		  define fileversion 1
		  define filetype $NPTYPE
                 }
		 if((temptemp=='cdum')||(temptemp=='icdu')){\
		  define fileversion 1
		  define filetype $CALCTYPE
                 }
		 if((temptemp=='fldum')||(temptemp=='ifldu')){\
		  define fileversion 1
		  define filetype $FIELDLINETYPE
                 }
		 if(('$1'=='0_avg2d.dat')||('$1'=='i0_avg2d.dat')){\
 		  define fileversion 1
		  define filetype $AVG2DTYPE
                 }
		 if('$1'=='0_avg1d.dat'){\
 		  define fileversion 1
		  define filetype $AVG1DTYPE
                 }
		 if('$1'=='0_ener.dat'){\
		  # ver=5  back to =n11 (inj/rad new)(n11 has no labels)
                  # ver=4  back to =n06 (inj new)
                  # ver=3  back to =o01
		  define fileversion 5
		  define filetype $ENERTYPE
                 }
		 if('$1'=='0_mode.dat'){\
		  define fileversion 1
		  define filetype $MODETYPE
                 }
		 if('$1'=='0_loss.dat'){\
                  # ver=5 back to o01
 		  define fileversion 5
		  define filetype $LOSSTYPE
                 }
		 if('$1'=='0_logsp.dat'){\
 		  define fileversion 1
		  define filetype $SPTYPE
                 }
		 if('$1'=='0_timescales.dat'){\
 		  define fileversion 1
		  define filetype $TSTYPE
                 }
		 if( ('$1'=='0_logdt.out.00')||('$1'=='0_logdt.out.01')||('$1'=='0_logdt.out.02')||('$1'=='0_logdt.out.03') ){\
 		  define fileversion 1
		  define filetype $LOGDTTYPE
                 }
		 if('$1'=='0_logstep.out'){\
 		  define fileversion 1
		  define filetype $STEPTYPE
                 }
		 if('$1'=='0_logperf.out'){\
		  define fileversion 3 # 2 is real, but 3 for now
		  define filetype $PERFTYPE
                 }
		 if((temptemp=='hst0')||(temptemp=='hst1')){\
 		  define fileversion 100
		  define filetype $ENERTYPE
                 }
                 define LSTART -1 define LFINAL 0
		}\
		else {\
                 gpa 0
                 define LSTART (1)
                 define LFINAL (2)
                 lines $LSTART $LFINAL
                 read {_fileversion 1 _filetype 2}
                 define fileversion (_fileversion)
                 define filetype (_filetype)
                }

                #
readpar     3   #
		readpar1 $1 $2 $3
		readpar11 $1 $2 $3
		readpar2 $1 $2 $3
		readpar3 $1 $2 $3
		#
readpar1     3  #
                # rev=1 starts with o1,o7
                # rev start with o01  and goes to o11 (no alpha and no averages and no mass)
                if($fileversion==1){\
		 # rev probably defined by alphareal0, but no fallback
                 set _alphareal0=0,0,1
                 set _alphareal0[0]=-100
                 mylines
                 read {_nx 1 _ny 2 _nz 3}
                 mylines
		 read {_Sx 1 _Sy 2 _Sz 3 _Lx 4 _Ly 5 _Lz 6}
                 mylines
		 read {_rg 1 _cour 2 _cour2 3 _css 4 _gam 5 _wgam 6 _dt 7 _t 8 _tf 9}
                 mylines
                 read {_nuvnr 1 _nul 2 _nuten 3 _alphareal0 4 _nreal 5}
                 mylines
                 read {_GRAVC 1}
                 mylines
		 read {_DTl 1 _DTd 2 _DTi 3 _DTloss 4 _DTfloor 5 _DTtimestep 6}
                 mylines
                 read {_res 1 _nush 2}
                 mylines
                 read {_vgx 1 _vgy 2 _vgz 3}
                 define parreadtest (_alphareal0[0])
		}
                # rev=2 starts with n01,o12, includes n,upto f12
                if($fileversion==2){\
                 # rev defined by extra var on rg line
                 set _numavg=0,0,1
                 set _numavg[0]=-100
                 mylines
                 read {_nx 1 _ny 2 _nz 3}
                 mylines
		 read {_Sx 1 _Sy 2 _Sz 3 _Lx 4 _Ly 5 _Lz 6}
                 mylines
                 read {_rg 1 _cour 2 _cour2 3 _css 4 _alpha 5 _gam 6 _dt 7 _tstart 8 _tf 9 _tavgi 10 _tavgf 11 _numavg 12}
                 mylines
                 read {_nuvnr 1 _nul 2 _nuten 3 _alphareal0 4 _nreal 5}
                 mylines
                 read {_GRAVC 1 _MASS 2}
                 mylines
		 read {_DTl 1 _DTd 2 _DTi 3 _DTloss 4 _DTfloor 5 _DTtimestep 6 _DTpd 7}
                 mylines
                 read {_res 1 _nush 2}
                 mylines
                 read {_vgx 1 _vgy 2 _vgz 3}
                 define parreadtest (_numavg[0])
		}
                # rev=3 starts with n27 (coolfact, alpha, etc. move around)
                if($fileversion==3){\
                 # rev defined by DTener
                 set _DTener=0,0,1
                 set _DTener[0]=-100
                 mylines
                 read {_nx 1 _ny 2 _nz 3}
                 mylines
		 read {_Sx 1 _Sy 2 _Sz 3 _Lx 4 _Ly 5 _Lz 6}
                 mylines
                 read {_rg 1 _cour 2 _cour2 3 _css 4 coolfact 5 _gam 6 _dt 7 _tstart 8 _tf 9 _tavgi 10 _tavgf 11 _numavg 12}
                 mylines
                 read {_nuvnr 1 _nul 2 _nuten 3 _alphareal0 4 _nreal 5}
                 mylines
                 read {_GRAVC 1 _MASS 2}
                 mylines
                 read {_DTl 1 _DTd 2 _DTi 3 _DTloss 4 _DTfloor 5 _DTtimestep 6 _DTpd 7 _DTener 8}
                 mylines
                 read {_res 1 _nush 2}
                 mylines
                 read {_vgx 1 _vgy 2 _vgz 3}
                 define parreadtest (_DTener[0])
		}
                # ver=4 starts with run f12, includes g
                if($fileversion==4){\
                 # rev defined by extra var on rg line
                 set _numavg=0,0,1
                 set _numavg[0]=-100
                 mylines
                 read {_nx 1 _ny 2 _nz 3}
                 mylines
		 read {_Sx 1 _Sy 2 _Sz 3 _Lx 4 _Ly 5 _Lz 6}
                 mylines
                 read {_rg 1 _rgp 2 _cour 3 _cour2 4 _css 5 coolfact 6 _gam 7 _dt 8 _tstart 9 _tf 10 _tavgi 11 _tavgf 12 _numavg 13}
                 mylines
                 read {_nuvnr 1 _nul 2 _nuten 3 _alphareal0 4 _nreal 5}
                 mylines
                 read {_GRAVC 1 _MASS 2}
                 mylines
                 read {_DTl 1 _DTd 2 _DTi 3 _DTloss 4 _DTfloor 5 _DTtimestep 6 _DTpd 7 _DTener 8}
                 mylines
                 read {_res 1 _nush 2}
                 mylines
                 read {_vgx 1 _vgy 2 _vgz 3}
                 define parreadtest (_numavg[0])
		}
                # ver=5 starts with a1
                # 5 is version as of Nov 1, 2000
                if($fileversion==5){\
                 # rev defined by extra line
                 set _vgz=0,0,1
                 set _vgz[0]=-100
                 mylines
                 read {_nx 1 _ny 2 _nz 3}
                 mylines
		 read {_Sx 1 _Sy 2 _Sz 3 _Lx 4 _Ly 5 _Lz 6}
                 mylines
                 read {_rg 1 _rgp 2 _css 3 coolfact 4 _gam 5 _alphareal0 6 _nreal 7}
                 mylines
                 read {_nuvnr 1 _nul 2 _nuten 3 _cour 4 _cour2 5}
                 mylines
                 read {_GRAVC 1 _MASS 2}
                 mylines
		 read {_tstart 1 _tf 2 _tavgi 3 _tavgf 4 _numavg 5}
                 mylines
                 read {_DTl 1 _DTd 2 _DTi 3 _DTloss 4 _DTfloor 5 _DTtimestep 6 _DTpd 7 _DTener 8}
                 mylines
                 read {_res 1 _nush 2}
                 mylines
                 read {_vgx 1 _vgy 2 _vgz 3}
                 define parreadtest (_vgz[0])
		}
		#
readpar11   3   #
                #
                # ver=6 starts with a8
                # 6 is version as of Nov 11, 2000
                if($fileversion==6){\
                 # rev defined by extra line
                 set _vgz=0,0,1
                 set _vgz[0]=-100
                 mylines
                 read {_nx 1 _ny 2 _nz 3}
                 mylines
		 read {_Sx 1 _Sy 2 _Sz 3 _Lx 4 _Ly 5 _Lz 6}
                 mylines
                 read {_rg 1 _rgp 2 _css 3 coolfact 4 _gam 5 _alphareal0 6 _nreal 7}
                 mylines
                 read {_nuvnr 1 _nul 2 _nuten 3 _cour 4 _cour2 5}
                 mylines
                 read {_GRAVC 1 _MASS 2}
                 mylines
		 read {_tstart 1 _tf 2 _tavgi 3 _tavgf 4 _numavg 5 _timagescale 6}
                 mylines
                 read {_DTl 1 _DTd 2 _DTi 3 _DTloss 4 _DTfloor 5 _DTtimestep 6 _DTpd 7 _DTener 8}
                 mylines
                 read {_res 1 _nush 2}
                 mylines
                 read {_vgx 1 _vgy 2 _vgz 3}
                 define parreadtest (_vgz[0])
		}
                #
                # ver=7 starts with mag field stuff
                # 7 is version as of Dec 10, 2000
                if($fileversion==7){\
                 # rev defined by extra lines
                 set _coord=0,0,1
                 set _coord[0]=-100
                 mylines
                 read {_nx 1 _ny 2 _nz 3}
                 mylines
		 read {_Sx 1 _Sy 2 _Sz 3 _Lx 4 _Ly 5 _Lz 6}
                 mylines
                 read {_rg 1 _rgp 2 _css 3 coolfact 4 _gam 5 _alphareal0 6 _nreal 7}
                 mylines
                 read {_nuvnr 1 _nul 2 _nuten 3 _cour 4 _cour2 5}
                 mylines
                 read {_GRAVC 1 _MASS 2}
                 mylines
		 read {_tstart 1 _tf 2 _tavgi 3 _tavgf 4 _numavg 5 _timagescale 6}
                 mylines
                 read {_DTl 1 _DTd 2 _DTi 3 _DTloss 4 _DTfloor 5 _DTtimestep 6 _DTpd 7 _DTener 8}
                 mylines
                 read {_res 1 _nush 2}
                 mylines
                 read {_vgx 1 _vgy 2 _vgz 3}
                 mylines
                 read {_coord 1 _fullvec 2 _analoutput 3}
                 mylines
                 read {_trans 1 _transx1 2 _transrhox1 3 _transiex1 4 _transv1x1 5 _transv2x1 6 _transv3x1 7 _transmagx1 8 _transx2 9 _transrhox2 10 _transiex2 11 _transv1x2 12 _transv2x2 13 _transv3x2 14}
                 mylines
                 read {_press 1 _pressx1 2 _pressx2 3}
                 mylines
                 read {_mag 1 _transbzx 2 _transbzy 3 _stepmocct 4 _mocctvx1 5 _mocctvx2 6 _mocctbx1 7 _mocctbx2 8}
                 mylines
                 read {_iedo 1 _viscart 2 _viscreal 3 _vreal 4 _vischeat 5 _mdotin 6 _cool 7 _res 8 _advint 9 _kever 10}
                 mylines
                 read {_intix1 1 _intox1 2 _intix2 3 _intox2 4 _intix3 5 _intox3 6}
                 mylines
                 read {_nonunigridx1 1 _nonunigridx2 2 _nonunigridx3 3 _simplebc 4 _bcix1 5 _bcox1 6 _bcix2 7 _bcox2 8 _bcix3 9 _bcox3 10}
                 #
                 define parreadtest (_coord[0])
		}
readpar2     3  #
                #
                # ver=8 starts after mhdtori colloqium fli
                # 8 is version as of Feb 17, 2001
                if($fileversion==8){\
                 # rev defined by extra lines
                 set _coord=0,0,1
                 set _coord[0]=-100
                 mylines
                 read {_nx 1 _ny 2 _nz 3}
                 mylines
		 read {_x1in 1 _x2in 2 _x3in 3 _x1out 4 _x2out 5 _x3out 6}
                 mylines
                 read {_rg 1 _rgp 2 _css 3 coolfact 4 _gam 5 _alphareal0 6 _nreal 7}
                 mylines
                 read {_nuvnr 1 _nul 2 _nuten 3 _cour 4 _cour2 5}
                 mylines
                 read {_GRAVC 1 _MASS 2}
                 mylines
		 read {_tstart 1 _tf 2 _tavgi 3 _tavgf 4 _numavg 5 _timagescale 6}
                 mylines
                 read {_DTl 1 _DTd 2 _DTi 3 _DTloss 4 _DTfloor 5 _DTtimestep 6 _DTpd 7 _DTener 8}
                 mylines
                 read {_res 1 _nush 2}
                 mylines
                 read {_vgx 1 _vgy 2 _vgz 3}
                 mylines
                 read {_coord 1 _fullvec 2 _analoutput 3}
                 mylines
                 read {_trans 1 _transx1 2 _transrhox1 3 _transiex1 4 _transv1x1 5 _transv2x1 6 _transv3x1 7 _transmagx1 8 _transx2 9 _transrhox2 10 _transiex2 11 _transv1x2 12 _transv2x2 13 _transv3x2 14}
                 mylines
                 read {_press 1 _pressx1 2 _pressx2 3}
                 mylines
                 read {_mag 1 _transbzx 2 _transbzy 3 _stepmocct 4 _mocctvx1 5 _mocctvx2 6 _mocctbx1 7 _mocctbx2 8}
                 mylines
                 read {_iedo 1 _viscart 2 _viscreal 3 _vreal 4 _vischeat 5 _mdotin 6 _cool 7 _res 8 _advint 9 _kever 10}
                 mylines
                 read {_intix1 1 _intox1 2 _intix2 3 _intox2 4 _intix3 5 _intox3 6}
                 mylines
                 read {_nonunigridx1 1 _nonunigridx2 2 _nonunigridx3 3 _simplebc 4 _bcix1 5 _bcox1 6 _bcix2 7 _bcox2 8 _bcix3 9 _bcox3 10}
                 #
                 set _Sx=_x1in
                 set _Sy=_x2in
                 set _Sz=_x3in
                 set _Lx=(_x1out-_x1in)
                 set _Ly=(_x2out-_x2in)
                 set _Lz=(_x3out-_x3in)
                 define parreadtest (_coord[0])
		}
                # ver=9 starts after gravitomagnetic tori
                # 9 is version as of March 12, 2001
                if($fileversion==9){\
                 # rev defined by extra lines
                 set _coord=0,0,1
                 set _coord[0]=-100
                 mylines
                 read {_nx 1 _ny 2 _nz 3}
                 mylines
		 read {_x1in 1 _x2in 2 _x3in 3 _x1out 4 _x2out 5 _x3out 6}
                 mylines
                 read {_rg 1 _rgp 2 _css 3 coolfact 4 _gam 5 _alphareal0 6 _nreal 7}
                 mylines
                 read {_nuvnr 1 _nul 2 _nuten 3 _cour 4 _cour2 5}
                 mylines
                 read {_GRAVC 1 _MASS 2 _invsol2 3 _blackholejz 4}
                 mylines
		 read {_tstart 1 _tf 2 _tavgi 3 _tavgf 4 _numavg 5 _timagescale 6}
                 mylines
                 read {_DTl 1 _DTd 2 _DTi 3 _DTloss 4 _DTfloor 5 _DTtimestep 6 _DTpd 7 _DTener 8}
                 mylines
                 read {_res 1 _nush 2}
                 mylines
                 read {_vgx 1 _vgy 2 _vgz 3}
                 mylines
                 read {_coord 1 _fullvec 2 _analoutput 3 _DYNAMICMM 4}
                 mylines
                 read {_trans 1 _transx1 2 _transrhox1 3 _transiex1 4 _transv1x1 5 _transv2x1 6 _transv3x1 7 _transmagx1 8 _transx2 9 _transrhox2 10 _transiex2 11 _transv1x2 12 _transv2x2 13 _transv3x2 14 _transvmagx3 15}
                 mylines
                 read {_press 1 _pressx1 2 _pressx2 3}
                 mylines
                 read {_mag 1 _transbzx 2 _transbzy 3 _stepmocct 4 _mocctvx1 5 _mocctvx2 6 _mocctbx1 7 _mocctbx2 8}
                 mylines
                 read {_iedo 1 _viscart 2 _viscreal 3 _vreal 4 _vischeat 5 _mdotin 6 _cool 7 _resreal 8 _rreal 9 _resheat 10 _advint 11 _kever 12}
                 mylines
                 read {_intix1 1 _intox1 2 _intix2 3 _intox2 4 _intix3 5 _intox3 6}
                 mylines
                 read {_nonunigridx1 1 _nonunigridx2 2 _nonunigridx3 3 _simplebc 4 _bcix1 5 _bcox1 6 _bcix2 7 _bcox2 8 _bcix3 9 _bcox3 10}
                 #
                 set _Sx=_x1in
                 set _Sy=_x2in
                 set _Sz=_x3in
                 set _Lx=(_x1out-_x1in)
                 set _Ly=(_x2out-_x2in)
                 set _Lz=(_x3out-_x3in)
                 define parreadtest (_coord[0])
		}
		#
readpar3 3      #
		# ver=10 starts after a bit of 3D implementation(before physics)
                # 10 is version as of April 3, 2001
                if($fileversion==10){\
                 # rev defined by extra lines
                 set _coord=0,0,1
                 set _coord[0]=-100
                 mylines
                 read {_nx 1 _ny 2 _nz 3}
                 mylines
		 read {_x1in 1 _x2in 2 _x3in 3 _x1out 4 _x2out 5 _x3out 6}
                 mylines
                 read {_rg 1 _rgp 2 _css 3 coolfact 4 _gam 5 _alphareal0 6 _nreal 7}
                 mylines
                 read {_nuvnr 1 _nul 2 _nuten 3 _cour 4 _cour2 5}
                 mylines
                 read {_GRAVC 1 _MASS 2 _invsol2 3 _blackholejz 4}
                 mylines
		 read {_tstart 1 _tf 2 _tavgi 3 _tavgf 4 _numavg 5 _timagescale 6}
                 mylines
                 read {_DTl 1 _DTd 2 _DTi 3 _DTloss 4 _DTfloor 5 _DTtimestep 6 _DTpd 7 _DTener 8 _DTfld 9}
                 mylines
                 read {_res 1 _nush 2}
                 mylines
                 read {_vgx 1 _vgy 2 _vgz 3}
                 mylines
                 read {_coord 1 _fullvec 2 _analoutput 3 _DYNAMICMM 4}
                 mylines
                 read {_trans 1 _transx1 2 _transrhox1 3 _transiex1 4 _transv1x1 5 _transv2x1 6 _transv3x1 7 _transmagx1 8 _transx2 9 _transrhox2 10 _transiex2 11 _transv1x2 12 _transv2x2 13 _transv3x2 14 _transmagx2 15 _transx3 16 _transrhox3 17 _transiex3 18 _transv1x3 19 _transv2x3 20 _transv3x3 21}
                 mylines
                 read {_press 1 _pressx1 2 _pressx2 3 _pressx3 4}
                 mylines
                 read {_mag 1 _transbzx 2 _transbzy 3 _stepmocct 4 _mocctvx1 5 _mocctvx2 6 _mocctvx3 7 _mocctbx1 8 _mocctbx2 9 _mocctbx3 10}
                 mylines
                 read {_iedo 1 _viscart 2 _viscreal 3 _vreal 4 _vischeat 5 _mdotin 6 _cool 7 _resreal 8 _rreal 9 _resheat 10 _advint 11 _kever 12}
                 mylines
                 read {_intix1 1 _intox1 2 _intix2 3 _intox2 4 _intix3 5 _intox3 6}
                 mylines
                 read {_nonunigridx1 1 _nonunigridx2 2 _nonunigridx3 3 _simplebc 4 _bcix1 5 _bcox1 6 _bcix2 7 _bcox2 8 _bcix3 9 _bcox3 10}
                 #
                 set _Sx=_x1in
                 set _Sy=_x2in
                 set _Sz=_x3in
                 set _Lx=(_x1out-_x1in)
                 set _Ly=(_x2out-_x2in)
                 set _Lz=(_x3out-_x3in)
                 define parreadtest (_coord[0])
		}
		#
		# ver=11 starts with DTmode
                # 11 is version as of July 13, 2001
                if($fileversion==11){\
                 # rev defined by extra lines
                 set _coord=0,0,1
                 set _coord[0]=-100
                 mylines
                 read {_nx 1 _ny 2 _nz 3}
                 mylines
		 read {_x1in 1 _x2in 2 _x3in 3 _x1out 4 _x2out 5 _x3out 6}
                 mylines
                 read {_rg 1 _rgp 2 _css 3 coolfact 4 _gam 5 _alphareal0 6 _nreal 7}
                 mylines
                 read {_nuvnr 1 _nul 2 _nuten 3 _cour 4 _cour2 5}
                 mylines
                 read {_GRAVC 1 _MASS 2 _invsol2 3 _blackholejz 4}
                 mylines
		 read {_tstart 1 _tf 2 _tavgi 3 _tavgf 4 _numavg 5 _timagescale 6}
                 mylines
                 read {_DTl 1 _DTd 2 _DTi 3 _DTloss 4 _DTfloor 5 _DTtimestep 6 _DTpd 7 _DTener 8 _DTfld 9 _DTmode 10}
                 mylines
                 read {_res 1 _nush 2}
                 mylines
                 read {_vgx 1 _vgy 2 _vgz 3}
                 mylines
                 read {_coord 1 _fullvec 2 _analoutput 3 _DYNAMICMM 4}
                 mylines
                 read {_trans 1 _transx1 2 _transrhox1 3 _transiex1 4 _transv1x1 5 _transv2x1 6 _transv3x1 7 _transmagx1 8 _transx2 9 _transrhox2 10 _transiex2 11 _transv1x2 12 _transv2x2 13 _transv3x2 14 _transmagx2 15 _transx3 16 _transrhox3 17 _transiex3 18 _transv1x3 19 _transv2x3 20 _transv3x3 21}
                 mylines
                 read {_press 1 _pressx1 2 _pressx2 3 _pressx3 4}
                 mylines
                 read {_mag 1 _transbzx 2 _transbzy 3 _stepmocct 4 _mocctvx1 5 _mocctvx2 6 _mocctvx3 7 _mocctbx1 8 _mocctbx2 9 _mocctbx3 10}
                 mylines
                 read {_iedo 1 _viscart 2 _viscreal 3 _vreal 4 _vischeat 5 _mdotin 6 _cool 7 _resreal 8 _rreal 9 _resheat 10 _advint 11 _kever 12}
                 mylines
                 read {_intix1 1 _intox1 2 _intix2 3 _intox2 4 _intix3 5 _intox3 6}
                 mylines
                 read {_nonunigridx1 1 _nonunigridx2 2 _nonunigridx3 3 _simplebc 4 _bcix1 5 _bcox1 6 _bcix2 7 _bcox2 8 _bcix3 9 _bcox3 10}
                 #
                 set _Sx=_x1in
                 set _Sy=_x2in
                 set _Sz=_x3in
                 set _Lx=(_x1out-_x1in)
                 set _Ly=(_x2out-_x2in)
                 set _Lz=(_x3out-_x3in)
                 define parreadtest (_coord[0])
		}
		#
                #
setvisc    0    #
                define fullvec 0
                define npdone 1
                set mag=0
                set _mag=0
                #
		#
rdbasic 3       #
		rdbasic1 $1 $2 $3 # parameters and setup
		rdbasic2 $1 $2 $3 # more setup
		rdbasic3 $1 $2 $3 # the grid itself
		rdbasic4 $1 $2 $3 # more setup
rdbasic1 3      #
                # [0/1][0/1]
                # [0][]->non interp [1][]->interp
                # [][0]->active grid [][1]->total grid
                # [][][n] number of cpu, or -1 for no multiple cpu
                # assume by default that does not exist unless specified after having done pp
                if($?npdone == 0){\
                  define npdone (0)
                }
                define rdbasicloaded (1)
                # number of files to skip for animation sequences
                define ANIMSKIP (1)
                set it=sprintf('%04d',$3)
                define cpunum (it)
                if($1==0){\
                 define filepar "0_gparam.par"
		 define interp (0)
		}\
		else{\
                 define filepar "0_igparam.par"
		 define interp (1)
		}
                if($2==0){\
                 define totalgrid (0)
		}\
		else{\
		 define totalgrid (1)
		}
                filetest $filepar
                readpar $1 $2 $3
                # check on 0_gparam.par read and redo if not good
                if($?parreadtest==1){\
                 if($parreadtest>-99){\
                  define temptemptemp (0)
                 }\
                 else{\
                  filetest $filepar
                  echo trying par fileversion 5
                  define fileversion 5
                  readpar $1 $2 $3
                 }
                 if($parreadtest>-99){\
                  define temptemptemp (0)
                 }\
                 else{\
                  filetest $filepar
                  echo trying par fileversion 4
                  define fileversion 4
                  readpar $1 $2 $3
                 }
                 if($parreadtest>-99){\
                  define temptemptemp (0)
                 }\
                 else{\
                  filetest $filepar
                  echo trying par fileversion 3
                  define fileversion 3
                  readpar $1 $2 $3
                 }
                 if($parreadtest>-99){\
                  define temptemptemp (0)
                 }\
                 else{\
                  filetest $filepar
                  echo trying par fileversion 2
                  define fileversion 2
                  readpar $1 $2 $3
                 }
                 if($parreadtest>-99){\
                  define temptemptemp (0)
                 }\
                 else{\
                  filetest $filepar
                  echo trying par fileversion 1
                  define fileversion 1
                  readpar $1 $2 $3
                 }
                }
		#
rdbasic2 3      #		
                #
                # define only needed vars
                define csgam1 (_css)
                define alphareal0 (_alphareal0)
		if($3>-1){\
                if($2==1){\
		       if(_nx>1){ define nx (_nx/$ncpux1+4) } else { define nx (_nx/$ncpux1) }
		       if(_ny>1){ define ny (_ny/$ncpux2+4) } else { define ny (_ny/$ncpux2) }
		       if(_nz>1){ define nz (_nz/$ncpux3+4) } else { define nz (_nz/$ncpux3) }
	        }\
                else{\
		       define nx (_nx/$ncpux1)
		       define ny (_ny/$ncpux2)
		       define nz (_nz/$ncpux3)
		       }
		    }
		if($3==-1){\
		if($2==1){\
		       if(_nx>1){ define nx (_nx+4) } else { define nx (_nx) }
		       if(_ny>1){ define ny (_ny+4) } else { define ny (_ny) }
		       if(_nz>1){ define nz (_nz+4) } else { define nz (_nz) }
	        }\
                else{\
		       define nx (_nx)
		       define ny (_ny)
		       define nz (_nz)
		       }
                }
		define gam (_gam)
                define fullvec (_fullvec)
                define coord (_coord)
                if(ABS($gam-1.0)<1.E-6){\
					  set wgam=0
					  set wgam1=1
					  set cs=1,$nx*$ny
					  set cs=$csgam1*cs/cs
		}\
		else{\
		       set wgam=1
		       set wgam1=0
		}
                if(ABS($gam-5.0/3.0)<1.E-6){\
					 set wgam53=1
		}\
		else{\
		       set wgam53=0
		}
                #
                #
rdbasic3 3      #		
                # now get geometry on active grid
                if($interp==0){\
                         if($totalgrid==0){\
				     if($cpunum<0){\
                                  rd 0_gridact1.par
                                  rd 0_gridact2.par 2
				     }\
				     else{\
                                  rd 0_gridact1.par.$cpunum
                                  rd 0_gridact2.par.$cpunum 2
					    }\
                         }\
                         else{\
				if($cpunum<0){\
					   rd 0_grid1.par
					   rd 0_grid2.par 2
				 }\
				    else{\
					   rd 0_grid1.par.$cpunum
					   rd 0_grid2.par.$cpunum 2
					   }\
			 }\
		}\
		else{\
                 if($totalgrid==0){\
				  if($cpunum<0){\
                                  rd 0_igridact1.par
                                  rd 0_igridact2.par 2
				     }\
				     else{\
                                  rd 0_igridact1.par.$cpunum
                                  rd 0_igridact2.par.$cpunum 2
					    }\
                         }\
                         else{\
				if($cpunum<0){\
					   rd 0_igridact1.par
					   rd 0_igridact2.par 2
				 }\
				    else{\
					   rd 0_igridact1.par.$cpunum
					   rd 0_igridact2.par.$cpunum 2
					   }\
			 }\
		}		
                #
rdbasic4 3	#	
                # define derivative vars
                # for full grid require data that global param file does not have
                if($2==1){\
		 if($nx>1){ define Lx (x12[$nx-1]-x12[0]) }\
		 else{ define Lx (x12[0])}
		 if($nx>1){ define Ly (x22[$nx*$ny-1]-x22[0]) }\
		 else{ define Ly (x22[0])}
		 if($nx>1){ define Lz (x3[$nx*$ny*$nz-1]-x3[0]) }\
		 else{ define Lz (x32[0])}
                 # 
		 define Sx (x12[0])
		 define Sy (x22[0]) 
		 define Sz (x32[0]) # 3D
		 #define Sz (x3[0])      # 2D
                }\
                else{\
                define Lx (_Lx) 
                define Ly (_Ly)
                define Lz (_Lz)
                define Sx (_Sx)
                define Sy (_Sy)
                define Sz (_Sz)
		       }
                define dx ($Lx/$nx)
		define dy ($Ly/$ny)
		define dz ($Lz/$nz)
                #
rdnumd 0        #
                da 0_numdumps.dat
                lines 1 2
                read {_NUMDUMPS 1}
                define NUMDUMPS (_NUMDUMPS)
		da 0_numfldumps.dat
                lines 1 2
                read {_NUMFLDUMPS 1}
                define NUMFLDUMPS (_NUMFLDUMPS)
                da 0_numpdumps.dat
                lines 1 2
                read {_NUMPDUMPS 1}
                define NUMPDUMPS (_NUMPDUMPS)
                da 0_numadumps.dat
                lines 1 2
                read {_NUMADUMPS 1}
                define NUMADUMPS (_NUMADUMPS)
                da 0_numnpdumps.dat
                lines 1 2
                read {_NUMNPDUMPS 1}
                define NUMNPDUMPS (_NUMNPDUMPS)
                da 0_numfloordumps.dat
                lines 1 2
                read {_NUMFLOORDUMPS 1}
                define NUMFLOORDUMPS (_NUMFLOORDUMPS)
		da ./i/0_numimages.dat
                lines 1 2
                read {_NUMIMAGES 1}
                define NUMIMAGES (_NUMIMAGES)
                echo filenumber=(0..NUM-1)
                echo Number of dumps: $NUMDUMPS
                echo Number of npdumps: $NUMNPDUMPS
                echo Number of pdumps: $NUMPDUMPS
                echo Number of adumps: $NUMADUMPS
                echo Number of floordumps: $NUMFLOORDUMPS
		echo Number of images: $NUMIMAGES
		echo Number of fldumps: $NUMFLDUMPS
myrd   12       # (checks filename, loads if name, does not if name is 0)
		 set tempit=WHATIS($1)
		 if(tempit==0){\
		  define temptemptemp (1)
		 }
		 if(tempit!=0) {
		  if($GAMMIE==1) {
		   jrdp $1
		  }
	          if($GAMMIE==0) {
	           if($?2 == 0) {
	            rd $1
	           }
	           if($?2 == 1) {
	            rd $1 $2
	           }
                  }
	         }
                 echo "rd done"
                 #
rd	12	# read macro for all data
                # e.g. rd 0_loss.dat  e.g. rd adump0000.dat 2
                if($?finaldraft==0){ rdraft }
                define failuremode (0)
                #
                if($?2 == 1){\
		 if($2==-1){\
		  define temptemptemp (0)
                  # just do not change anything!
                  # needed since otherwise rd np$1 inside would switch definitions
                 }
		 if($2==-2){\
 		  define exttype (0)
                  define ext " "
                  define extx x
                  define exty y
                  define extz z
                 }
                 if(($2!=-1)&&($2!=-2)){\
		  define exttype (1)
		  define ext $2
		  define extx $2x
                  define exty $2y
                  define extz $2z
		 }
                }\
                else{\
		 define exttype (0)
                 define ext " "
                 define extx x
                 define exty y
                 define extz z
                }
                rdintro $1
		if(!($failuremode)){\
                 rdreads $1
                 if($filetype==$ENERTYPE){\
                  # check on 0_ener.dat read and redo if not good
                  if($?enerreadtest==1){\
                   if($enerreadtest>-99){\
                    define temptemptemp (0)
                   }\
                   else{\
                    echo trying ener fileversion 4
                    rdintro $1
                    define fileversion 4
                    rdreads $1
                   }
                  }
                  if($?enerreadtest==1){\
                   if($enerreadtest>-99){\
                    define temptemptemp (0)
                   }\
                   else{\
                    echo trying ener fileversion 3
                    rdintro $1
                    define fileversion 3
                    rdreads $1
                   }
                  }
                 }
                 if($filetype==$NPTYPE){\
 		  if($npdumpreadtest<-99){\
                   echo trying npdump fileversion 0
                   rdintro $1
                   define fileversion 0
		   rdreads $1
		  }
		 }
                 #
                 if( ($filetype==$DTYPE)&&($npdone) ){\
                  # generally read in npdump when dumps are read
	 	  define oldfileversion2 $fileversion
		  define oldfiletype2 $filetype
	 	  rd np$1 -1
                  # return to normal filetype after done
                  define fileversion $oldfileversion2
                  define filetype $oldfiletype2
		 }
                 # assume dump gets the calc for npdump
                 if($filetype!=$NPTYPE){\
                  if($DOCALCS){\
                   rdcalcs1 $1
                   rdcalcs2 $1
                  }
		 }
		}
                #
rdbondi   0     # (must load rdbasic still first for gam)(why the fuck is this if statement not same as others in back slashes!)
                da 0_analdata2.dat
                read {rad 1 rho 2 ie 3 gravpot 4 vr 5 vtheta 6 vphi 7 br 8 btheta 9 bphi 10}
                #
                if(wgam1==1){
		  set cs=0,($nx+4)*$ny-1
		  do ii=0,($nx+4)*$ny-1,1 {
		    set cs[$ii]=$csgam1
		  }
			       set cs2=cs*cs
			       set Be=0.5*vr*vr+cs*cs*LN(rho)+gravpot
		}\
		else{
		       set cs2=$gam*($gam-1)*ie/rho
                        set cs=sqrt(cs2)

		       set Be=0.5*vr*vr+cs*cs/($gam-1)+gravpot                  
		}
                set Machr=vr/cs
		set Mdot=4*3.14159*rad*rad*rho*vr
		set p=($gam-1.0)*ie
		# entropy per degree of freedom
		set entropy=p/rho**($gam)
		set me=9.1*10**(-28)
		set hbar=1.05459*10**(-27)
		set kb=1.38066*10**(-16)
		set electronentropy=entropy*(me**($gam+1))/(2*pi*kb*hbar**2)
                #
magcall         #
		set mflux=(r*vx* x1*x1*sin(x2)*dx2)
		set keflux=(r*vx*0.5*(vx**2+vy**2+vz**2)*x1*x1*sin(x2)*dx2)
                set hflux = (vx*($gam*en)*x1*x1*sin(x2)*dx2)
		set beflux=((bx**2+by**2+bz**2)*vx-(bx*vx+by*vy+bz*vz)*bx)*x1*x1*sin(x2)*dx2
		set peflux=(r*vx*pot*x1*x1*sin(x2)*dx2)
		set eflux = keflux+hflux+beflux+peflux
		set v1flux=(x1*x1*sin(x2)*dx2)*(r*vx*vx+p-bx*bx+0.5*(bx**2+by**2+bz**2))
		set v3fluxrey=(r*vx*vz* x1*x1*sin(x2)*dx2 * x1*sin(x2))
		set v3fluxmag=-(bx*bz* x1*x1*sin(x2)*dx2 * x1*sin(x2))
                set v3flux = v3fluxrey+v3fluxmag
		set b1flux2=x1*sin(x2)*(vx*by-vy*bx)*dx1
		set b2flux1=-x1*sin(x2)*(vx*by-vy*bx)*dx2
		set b3flux1=x1*(vx*bz-vz*bx)*dx2
		set b3flux2=-(vy*bz-vz*vy)*dx1
                # magbondi
                define mag (1) define npdone (0) define coord (1)
                #set MAGK = 0.3317559754612477089752102542397933695730641977564196
                #set MAGF = 8.4822847861354036832817652116918313741194889910534187E-10
                #set MAGOMEGA = 0.000044721359549995793928183473374625524708812367192210704
                #set MAGPHI = 0.0095086263442549441567086161486063491924189885097233
                #set xs = 0.7770330780030363194074308952074505
                #set ys = 1.940544269119426993709021901238330
                #set xf = 1.302391299892080014807550863286971
                #set yf = 0.5139371598659040285268997762736311
                # set En = 1.738429467866434807018809306133664 // dimensionless form
                #set MAGEN = 0.003476858935732869614037618612267327
                #set beta = 0.5755932415448330480108551469612932
                set MAGRA = 500.0
                set RHORA = 1E-13
                set gam=1.2
                set GM = 1
		set MAGTHETA=(0.5) # // dimenless
                set MAGomega=(.25) # // dimenless
		set  SPFAST=(1.302391299892080014807550863286971) # // y=0.5139371598659040285268997762736311
                set SPALF=(1.0) #// y=1.0
                set SPSLOW=(0.7770330780030363194074308952074505) #// y=1.940544269119426993709021901238330
                set MAGBETA=(0.5755932415448330480108551469612932) #// dimenless
                set MAGE=(1.738429467866434807018809306133664) #// dimenless:
		set MAGK=(RHORA**(1.0-gam)*MAGTHETA/(gam*MAGRA))
                set MAGF=(RHORA*SQRT(MAGBETA*GM*MAGRA**(3.0)))
                set MAGOMEGA=(SQRT(GM*MAGomega/MAGRA**(3.0)))
                # #define MAGPHI (sqrt(4.0*M_PI*MAGBETA*GM*pow(MAGRA,3.0)*RHORA))
                set MAGPHI=(SQRT(MAGBETA*GM*MAGRA**(3.0)*RHORA))
                set VELRA=(MAGF/(MAGRA*MAGRA*RHORA))
                set MAGENERGY=(MAGE*GM/MAGRA)
                # mag bondi equations for check, should all be 0
                rdbondi
                set p=($gam-1.0)*ie
                set diff1=p-MAGK*rho**($gam)
                set diff2=MAGF-rho*vr*rad*rad
                set diff3=MAGPHI-br*rad*rad
                set diff4=vr*bphi-br*(vphi-MAGOMEGA*rad)
                set diff5=MAGOMEGA*MAGRA*MAGRA-rad*(vphi-br*bphi/(rho*vr))
                set diff6=0.5*(vr**2+(vphi-MAGOMEGA*rad)**2)+$gam/($gam-1)*p/rho-GM/rad-0.5*MAGOMEGA*MAGOMEGA*rad*rad-MAGENERGY
                #set ftemp=rad*(MAGPHI*MAGPHI*rho-MAGF*MAGF)
                #set vphinum=MAGOMEGA*(MAGPHI*MAGPHI*rad*rad*rho-MAGF*MAGF*MAGRA*MAGRA)
                set energy=1/2*rho*(vr**2+vtheta**2+vphi**2)+ie+1/(2)*(br**2+btheta**2+bphi**2)+rho*gravpot
		set senergy=1/2*(vr**2+vtheta**2+vphi**2)+ie/rho+1/(2*rho)*(br**2+btheta**2+bphi**2)+rho*gravpot
                #
rdintro 1       #
		da $1
                filetest $1
                if(($filetype==$FIELDLINETYPE)||($filetype==$CALCTYPE)||($filetype==$DTYPE)||($filetype==$PDTYPE)||($filetype==$ADTYPE)||($filetype==$FLTYPE)||($filetype==$NPTYPE)||($filetype==$AVG2DTYPE)||($filetype==$AVG1DTYPE)){\
                 define LSTART ($LFINAL+1)
                 define LFINAL ($LSTART+1)
                 lines $LSTART $LFINAL
 	 	 read {_time 1 _SAMPLE 2 _ALLZONEC 3}
                 define time (_time)
                 define SAMPLE (_SAMPLE)
                 define ALLZONEC (_ALLZONEC)
                 define LSTART ($LFINAL+1)
                 define LFINAL ($LSTART+1)
		 lines $LSTART $LFINAL
  		 read {_N1 1 _N2 2 _N3 3}
                 define dnx (_N1)
                 define dny (_N2)
                 define dnz (_N3)
		 if($?rdbasicloaded==1){\
		  if( ($dnx!=$nx)||($dny!=$ny)||($dnz!=$nz) ){\
		   echo dump grid different than grid size expects, load rdbasic?
		   define failuremode (1)
		   return
		  }
		 }\
		 else{\
                  define oldLSTART $LSTART
                  define oldLFINAL $LFINAL
		  define oldexttype $exttype
		  define oldext "$!ext"
		  define oldextx "$!extx"
                  define oldexty "$!exty"
                  define oldextz "$!extz"
		  define oldfileversion $fileversion
		  define oldfiletype $filetype
                  rdbasic 0 0 -1 # load typical params
                  define LSTART $oldLSTART
                  define LFINAL $oldLFINAL
		  define exttype $oldexttype
		  define ext "$!oldext"
		  define extx "$!oldextx"
                  define exty "$!oldexty"
                  define extz "$!oldextz"
                  define fileversion $oldfileversion
                  define filetype $oldfiletype
		 }
                }
                #
rdreads 1       #
                rdreads1 $1
                rdreads2 $1
                rdreads22 $1
                rdreads3 $1
		rdreads31 $1
		rdreads32 $1
                rdreads4 $1
                #
rdreads1 1      #
                da $1   # in case other reads have occured so far
                # generally good
                define LSTART ($LFINAL+1)
                define LFINAL 1000000
		lines $LSTART $LFINAL
                #echo thelines $LSTART $LFINAL
                #con$!!ext 11
                if(($filetype==$DTYPE)||($filetype==$ADTYPE)||($filetype==$PDTYPE)){\
		 if($fullvec){\
		  read {r$!!ext 1 en$!!ext 2 pot$!!ext 3 v$!!extx 4 v$!!exty 5 v$!!extz 6 diag1$!!ext 7 diag2$!!ext 8 diag3$!!ext 9 diag4$!!ext 10 diag5$!!ext 11 diag6$!!ext 12 diag7$!!ext 13 diag8$!!ext 14 diag9$!!ext 15 diag10$!!ext 16 diag11$!!ext 17 diag12$!!ext 18 diag13$!!ext 19 diag14$!!ext 20 diag15$!!ext 21 diag16$!!ext 22 diag17$!!ext 23 diag18$!!ext 24 diag19$!!ext 25 diag20$!!ext 26 diag21$!!ext 27 diag22$!!ext 28 diag23$!!ext 29 diag24$!!ext 30 diag25$!!ext 31 diag26$!!ext 32 diag27$!!ext 33 diag28$!!ext 34 diag29$!!ext 35 diag30$!!ext 36 diag31$!!ext 37 diag32$!!ext 38 diag33$!!ext 39 }
		}\
		 else{\
		  read {r$!!ext 1 en$!!ext 2 pot$!!ext 3 v$!!extx 4 v$!!exty 5 v$!!extz 6}
		 }
                }
                if($filetype==$CALCTYPE){\
		 read {LR$!!extx 1 LR$!!exty 2 LR$!!extz 3 LV$!!extx 4 LV$!!exty 5 LV$!!extz 6 ER$!!extx 7 ER$!!exty 8 ER$!!extz 9 EV$!!extx 10 EV$!!exty 11 EV$!!extz 12}
                }
		if($filetype==$FIELDLINETYPE){\
		 read {phib$!!extz 1}
                }
                if($filetype==$FLTYPE){\
                 if($fileversion==2){\
   		  read {r$!!ext 1 en$!!ext 2 pot$!!ext 3 ke$!!ext 4 v$!!extx 5 v$!!exty 6 v$!!extz 7}
                 }
                 if($fileversion==1){\
   		  read {r$!!ext 1 en$!!ext 2 pot$!!ext 3 v$!!extx 4 v$!!exty 5 v$!!extz 6}
                 }
                }
                if($filetype==$NPTYPE){\
					 if(0){\
                 set nuvisc$!!ext=0,0,1
                 set nuvisc$!!ext[0]=-100
                 if($fileversion==2){\
		  read {sig11$!!ext 1 sig12$!!ext 2 sig13$!!ext 3 sig22$!!ext 4 sig23$!!ext 5 sig33$!!ext 6 nuvisc$!!ext 7}
                  set sig21$!!ext=sig12$!!ext
                  set sig31$!!ext=sig13$!!ext
                  set sig32$!!ext=sig23$!!ext
                 }
                 if($fileversion==1){\
		  read {sig11$!!ext 1 sig12$!!ext 2 sig13$!!ext 3 sig21$!!ext 4 sig22$!!ext 5 sig23$!!ext 6 sig31$!!ext 7 sig32$!!ext 8 sig33$!!ext 9 nuvisc$!!ext 10}
                 }
                 if($fileversion==0){\
 		  read {sig11$!!ext 1 sig12$!!ext 2 sig13$!!ext 3 sig21$!!ext 4 sig22$!!ext 5 sig23$!!ext 6 sig31$!!ext 7 sig32$!!ext 8 sig33$!!ext 9}
                 }
                 define npdumpreadtest (nuvisc$!!ext[0])
		   }
		if(1){\
			read {avx 1 avy 2 avz 3 apx 4 apy 5 apz 6 agx 7 agy 8 agz 9 acbx 10 acby 11 acbz 12 ab2x 13 ab2y 14 ab2z 15 abx 16 aby 17 abz 18 \
				fkex 19 fkey 20 fkez 21 fhx 22 fhy 23 fhz 24 fgx 25 fgy 26 fgz 27 fbx 28 fby 29 fbz 30 fb2x 31 fb2y 32 fb2z 33 fmx 34 fmy 35 fmz 36 fpxx 37 fpxy 38 fpxz 39 fpyx 40 fpyy 41 fpyz 42 fpzx 43 fpzy 44 fpzz 45 fbxx 46 fbxy 47 fbxz 48 fbyx 49 fbyy 50 fbyz 51 fbzx 52 fbzy 53 fbzz 54 \
				vm 55 vkex 56 vkey 57 vkez 58 vie 59 vge 60 vbex 61 vbey 62 vbez 63 vpx 64 vpy 65 vpz 66 vlx 67 vly 68 vlz 69 vbx 70 vby 71 vbz 72}
		set atotx=avx+apx+agx+ab2x+abx
		set atoty=avy+apy+agy+ab2y+aby
		set atotz=avz+apz+agz+ab2z+abz
		set amagx=ab2x+abx
		set amagy=ab2y+aby
		set amagz=ab2z+abz
		}
		define npdumpreadtest (avx[0])
		  }
                if($filetype==$AVG2DTYPE){\
                 if($fileversion==1){\
                  read {r$!!ext 1 en$!!ext 2 be2d$!!ext 3 csq2d$!!ext 4 e2d$!!ext 5 v$!!extx 6 v$!!exty 7 v$!!extz 8}
                 }
                 if($fileversion==2){\
		  read {r$!!ext 1 en$!!ext 2 be2d$!!ext 3 csq2d$!!ext 4 e2d$!!ext 5 v$!!extx 6 v$!!exty 7 v$!!extz 8 sig11$!!ext 9 sig12$!!ext 10 sig13$!!ext 11 sig22$!!ext 12 sig23$!!ext 13 sig33$!!ext 14 nuvisc$!!ext 15}
                  set sig21$!!ext=sig12$!!ext
                  set sig31$!!ext=sig13$!!ext
                  set sig32$!!ext=sig23$!!ext
                 }		
                }
                if($filetype==$AVG1DTYPE){\
                 if($fileversion==1){\
  		  read {r1d1$!!ext 1 en1d1$!!ext 2 be1d1$!!ext 3 csq1d1$!!ext 4 s1d1$!!ext 5 v1d1$!!extx 6 v1d1$!!exty 7 v1d1$!!extz 8\
		        r1d2$!!ext 9 en1d2$!!ext 10 be1d2$!!ext 11 csq1d2$!!ext 12 s1d2$!!ext 13 v1d2$!!extx 14 v1d2$!!exty 15 v1d2$!!extz 16}
                 }
                 if($fileversion==2){\
  	          read {r1d1$!!ext 1 en1d1$!!ext 2 be1d1$!!ext 3 csq1d1$!!ext 4 s1d1$!!ext 5 v1d1$!!extx 6 v1d1$!!exty 7 v1d1$!!extz 8 sig111d1$!!ext 9 sig121d1$!!ext 10 sig131d1$!!ext 11 sig221d1$!!ext 12 sig231d1$!!ext 13 sig331d1$!!ext 14 nuvisc1d1$!!ext 15\
	       	        r1d2$!!ext 16 en1d2$!!ext 17 be1d2$!!ext 18 csq1d2$!!ext 19 s1d2$!!ext 20 v1d2$!!extx 21 v1d2$!!exty 22 v1d2$!!extz 23 sig111d2$!!ext 24 sig121d2$!!ext 25 sig131d2$!!ext 26 sig221d2$!!ext 27 sig231d2$!!ext 28 sig331d2$!!ext 29 nuvisc1d2$!!ext 30}
	         }
                  echo make sure you read in normal dump for potential
  		  if($AVG1DWHICH==1){\
                   set r$!!ext=r1d1$!!ext
                   set en$!!ext=en1d1$!!ext
                   set v$!!extx=v1d1$!!extx
                   set v$!!exty=v1d1$!!exty
                   set v$!!extz=v1d1$!!extz
                   set sig11$!!ext=sig131d1$!!ext
                   set sig12$!!ext=sig121d1$!!ext
                   set sig13$!!ext=sig131d1$!!ext
                   set sig22$!!ext=sig221d1$!!ext
                   set sig23$!!ext=sig231d1$!!ext
                   set sig33$!!ext=sig331d1$!!ext
	          }
  	          if($AVG1DWHICH==2){\
                   set r$!!ext=r1d2$!!ext
                   set en$!!ext=en1d2$!!ext
                   set v$!!extx=v1d2$!!extx
                   set v$!!exty=v1d2$!!exty
                   set v$!!extz=v1d2$!!extz
                   set sig11$!!ext=sig131d2$!!ext
                   set sig12$!!ext=sig121d2$!!ext
                   set sig13$!!ext=sig131d2$!!ext
                   set sig22$!!ext=sig221d2$!!ext
                   set sig23$!!ext=sig231d2$!!ext
                   set sig33$!!ext=sig331d2$!!ext
	          }
                }
                if($filetype==$GRIDTYPE){\
                 if($fileversion==2){\
                   read {gr$!!ext 1 k$!!ext 2 j$!!ext 3 i$!!ext 4 dx1$!!ext 5 dx2$!!ext 6 dx3$!!ext 7 x1$!!ext 8 x2$!!ext 9 x3$!!ext 10 g1$!!ext 11 dg1$!!ext 12 g2$!!ext 13 dg2$!!ext 14\
	 	         g3$!!ext 15 dg3$!!ext 16 g4$!!ext 17 dg4$!!ext 18 cot$!!ext 19 dvl1$!!ext 20 dvl2$!!ext 21 dvl3$!!ext 22\
		         bcstyp1$!!ext 23 bcsdim1$!!ext 24 bcsdir1$!!ext 25 bcstyp2$!!ext 26 bcsdim2$!!ext 27 bcsdir2$!!ext 28 bcstyp3$!!ext 29 bcsdim3$!!ext 30 bcsdir3$!!ext 31\
		         bcvtyp1$!!ext 32 bcvdim1$!!ext 33 bcvdir1$!!ext 34 bcvtyp1$!!ext 35 bcvdim1$!!ext 36 bcvdir1$!!ext 37}
                   set xc$!!ext=x1$!!ext*sin(x2$!!ext)
                   set zc$!!ext=x1$!!ext*cos(x2$!!ext)
                 }
		 if($fileversion==3){\
                   read {gr$!!ext 1 k$!!ext 2 j$!!ext 3 i$!!ext 4 dx1$!!ext 5 dx2$!!ext 6 dx3$!!ext 7 x1$!!ext 8 x2$!!ext 9 x3$!!ext 10 g1$!!ext 11 dg1$!!ext 12 g2$!!ext 13 dg2$!!ext 14\
	 	         g3$!!ext 15 dg3$!!ext 16 g4$!!ext 17 dg4$!!ext 18 dvl1$!!ext 19 dvl2$!!ext 20 dvl3$!!ext 21\
		         bcstyp1$!!ext 22 bcsdim1$!!ext 23 bcsdir1$!!ext 24 bcstyp2$!!ext 25 bcsdim2$!!ext 26 bcsdir2$!!ext 27 bcstyp3$!!ext 28 bcsdim3$!!ext 29 bcsdir3$!!ext 30\
		         bcvtyp1$!!ext 31 bcvdim1$!!ext 32 bcvdir1$!!ext 33 bcvtyp1$!!ext 34 bcvdim1$!!ext 35 bcvdir1$!!ext 36}
                   set xc$!!ext=x1$!!ext*sin(x2$!!ext)
                   set zc$!!ext=x1$!!ext*cos(x2$!!ext)
                 }
                }
                #
rdreads2  1     #
                if($filetype==$ENERTYPE){\
                 if($fileversion==7){\
                  set mrad$!!ext=0,0,1
                  set mrad$!!ext[0]=-100
		  if($fullvec){\
		   read {t$!!ext 1\
                   etot$!!ext 2 etotdel$!!ext 3 etotblost$!!ext 4 etotfladd$!!ext 5 etotinj$!!ext 6 etotrad$!!ext  7 etotdiff$!!ext  8\
		   m$!!ext    9 mdel$!!ext   10 mblost$!!ext   11 mfladd$!!ext   12 minj$!!ext   13 mrad$!!ext    14 mdiff$!!ext    15\
                   ie$!!ext  16 iedel$!!ext  17 ieblost$!!ext  18 iefladd$!!ext  19 ieinj$!!ext  20 ierad$!!ext   21 iediff$!!ext   22\
                   pe$!!ext 23 pedel$!!ext 24 peblost$!!ext 25 pefladd$!!ext 26 peinj$!!ext 27 perad$!!ext  28 pediff$!!ext  29\
		   ke$!!ext  30 kedel$!!ext  31 keblost$!!ext  32 kefladd$!!ext  33 keinj$!!ext  34 kerad$!!ext   35 kediff$!!ext   36\
                   mx1$!!ext 37 mx1del$!!ext 38 mx1blost$!!ext 39 mx1fladd$!!ext 40 mx1inj$!!ext 41 mx1rad$!!ext  42 mx1diff$!!ext  43\
                   mx2$!!ext 44 mx2del$!!ext 45 mx2blost$!!ext 46 mx2fladd$!!ext 47 mx2inj$!!ext 48 mx2rad$!!ext  49 mx2diff$!!ext  50\
		   mx3$!!ext 51 mx3del$!!ext 52 mx3blost$!!ext 53 mx3fladd$!!ext 54 mx3inj$!!ext 55 mx3rad$!!ext  56 mx3diff$!!ext  57\
                   be$!!ext  58 bedel$!!ext  59 beblost$!!ext  60 befladd$!!ext  61 beinj$!!ext  62 berad$!!ext   63 bediff$!!ext   64\
                   bx1$!!ext 65 bx1del$!!ext 66 bx1blost$!!ext 67 bx1fladd$!!ext 68 bx1inj$!!ext 69 bx1rad$!!ext  70 bx1diff$!!ext  71\
                   bx2$!!ext 72 bx2del$!!ext 73 bx2blost$!!ext 74 bx2fladd$!!ext 75 bx2inj$!!ext 76 bx2rad$!!ext  77 bx2diff$!!ext  78\
		   bx3$!!ext 79 bx3del$!!ext 80 bx3blost$!!ext 81 bx3fladd$!!ext 82 bx3inj$!!ext 83 bx3rad$!!ext  84 bx3diff$!!ext  85\
		   ve$!!ext  86 vedel$!!ext  87 veblost$!!ext  88 vefladd$!!ext  89 veinj$!!ext  90 verad$!!ext   91 vediff$!!ext   92\
		   amx1$!!ext 93 amx1del$!!ext 94 amx1blost$!!ext 95 amx1fladd$!!ext 96 amx1inj$!!ext 97 amx1rad$!!ext  98 amx1diff$!!ext  99\
                   amx2$!!ext 100 amx2del$!!ext 101 amx2blost$!!ext 102 amx2fladd$!!ext 103 amx2inj$!!ext 104 amx2rad$!!ext  105 amx2diff$!!ext  106\
		   amx3$!!ext 107 amx3del$!!ext 108 amx3blost$!!ext 109 amx3fladd$!!ext 110 amx3inj$!!ext 111 amx3rad$!!ext  112 amx3diff$!!ext  113}
                   #
		  }
		  # no $fullvec==0 case -- too annoying to redo the numbers
		  define enerreadtest (etotdel$!!ext[0])
                 }
                 if($fileversion==6){\
                  set mrad$!!ext=0,0,1
                  set mrad$!!ext[0]=-100
		  if($fullvec){\
		   read {t$!!ext 1\
                   etot$!!ext 2 etotdel$!!ext 3 etotblost$!!ext 4 etotfladd$!!ext 5 etotinj$!!ext 6 etotrad$!!ext  7 etotdiff$!!ext  8\
		   m$!!ext    9 mdel$!!ext   10 mblost$!!ext   11 mfladd$!!ext   12 minj$!!ext   13 mrad$!!ext    14 mdiff$!!ext    15\
                   ie$!!ext  16 iedel$!!ext  17 ieblost$!!ext  18 iefladd$!!ext  19 ieinj$!!ext  20 ierad$!!ext   21 iediff$!!ext   22\
                   pe$!!ext 23 pedel$!!ext 24 peblost$!!ext 25 pefladd$!!ext 26 peinj$!!ext 27 perad$!!ext  28 pediff$!!ext  29\
		   ke$!!ext  30 kedel$!!ext  31 keblost$!!ext  32 kefladd$!!ext  33 keinj$!!ext  34 kerad$!!ext   35 kediff$!!ext   36\
                   mx1$!!ext 37 mx1del$!!ext 38 mx1blost$!!ext 39 mx1fladd$!!ext 40 mx1inj$!!ext 41 mx1rad$!!ext  42 mx1diff$!!ext  43\
                   mx2$!!ext 44 mx2del$!!ext 45 mx2blost$!!ext 46 mx2fladd$!!ext 47 mx2inj$!!ext 48 mx2rad$!!ext  49 mx2diff$!!ext  50\
		   mx3$!!ext 51 mx3del$!!ext 52 mx3blost$!!ext 53 mx3fladd$!!ext 54 mx3inj$!!ext 55 mx3rad$!!ext  56 mx3diff$!!ext  57\
                   be$!!ext  58 bedel$!!ext  59 beblost$!!ext  60 befladd$!!ext  61 beinj$!!ext  62 berad$!!ext   63 bediff$!!ext   64\
                   bx1$!!ext 65 bx1del$!!ext 66 bx1blost$!!ext 67 bx1fladd$!!ext 68 bx1inj$!!ext 69 bx1rad$!!ext  70 bx1diff$!!ext  71\
                   bx2$!!ext 72 bx2del$!!ext 73 bx2blost$!!ext 74 bx2fladd$!!ext 75 bx2inj$!!ext 76 bx2rad$!!ext  77 bx2diff$!!ext  78\
		   bx3$!!ext 79 bx3del$!!ext 80 bx3blost$!!ext 81 bx3fladd$!!ext 82 bx3inj$!!ext 83 bx3rad$!!ext  84 bx3diff$!!ext  85\
		   ve$!!ext  86 vedel$!!ext  87 veblost$!!ext  88 vefladd$!!ext  89 veinj$!!ext  90 verad$!!ext   91 vediff$!!ext   92}
                   #
		  }
		  # no $fullvec==0 case -- too annoying to redo the numbers
		  define enerreadtest (etotdel$!!ext[0])
                 }
		}
rdreads22    1  #
		if($filetype==$ENERTYPE){\
                 if($fileversion==5){\
                  set mrad$!!ext=0,0,1
                  set mrad$!!ext[0]=-100
		  if($fullvec){\
		   read {t$!!ext 1 cmode_amp$!!ext 2 smode_amp$!!ext 3\
                   etot$!!ext 4 etotdel$!!ext 5 etotblost$!!ext 6 etotfladd$!!ext 7 etotinj$!!ext 8 etotrad$!!ext 9 etotdiff$!!ext  10\
		   m$!!ext   11 mdel$!!ext   12 mblost$!!ext   13 mfladd$!!ext   14 minj$!!ext   15 mrad$!!ext    16 mdiff$!!ext    17\
                   ie$!!ext  18 iedel$!!ext  19 ieblost$!!ext  20 iefladd$!!ext  21 ieinj$!!ext  22 ierad$!!ext   23 iediff$!!ext   24\
                   pe$!!ext 25 pedel$!!ext 26 peblost$!!ext 27 pefladd$!!ext 28 peinj$!!ext 29 perad$!!ext  30 pediff$!!ext  31\
		   ke$!!ext  32 kedel$!!ext  33 keblost$!!ext  34 kefladd$!!ext  35 keinj$!!ext  36 kerad$!!ext   37 kediff$!!ext   38\
                   mx1$!!ext 39 mx1del$!!ext 40 mx1blost$!!ext 41 mx1fladd$!!ext 42 mx1inj$!!ext 43 mx1rad$!!ext  44 mx1diff$!!ext  45\
                   mx2$!!ext 46 mx2del$!!ext 47 mx2blost$!!ext 48 mx2fladd$!!ext 49 mx2inj$!!ext 50 mx2rad$!!ext  51 mx2diff$!!ext  52\
		   mx3$!!ext 53 mx3del$!!ext 54 mx3blost$!!ext 55 mx3fladd$!!ext 56 mx3inj$!!ext 57 mx3rad$!!ext  58 mx3diff$!!ext  59\
                   be$!!ext  60 bedel$!!ext  61 beblost$!!ext  62 befladd$!!ext  63 beinj$!!ext  64 berad$!!ext   65 bediff$!!ext   66\
                   bx1$!!ext 67 bx1del$!!ext 68 bx1blost$!!ext 69 bx1fladd$!!ext 70 bx1inj$!!ext 71 bx1rad$!!ext  72 bx1diff$!!ext  73\
                   bx2$!!ext 74 bx2del$!!ext 75 bx2blost$!!ext 76 bx2fladd$!!ext 77 bx2inj$!!ext 78 bx2rad$!!ext  79 bx2diff$!!ext  80\
		   bx3$!!ext 81 bx3del$!!ext 82 bx3blost$!!ext 83 bx3fladd$!!ext 84 bx3inj$!!ext 85 bx3rad$!!ext  86 bx3diff$!!ext  87\
		   ve$!!ext  88 vedel$!!ext  89 veblost$!!ext  90 vefladd$!!ext  91 veinj$!!ext  92 verad$!!ext   93 vediff$!!ext   94}
                   #
		  }\
		  else{\
		   read {t$!!ext 1 cmode_amp$!!ext 2 smode_amp$!!ext 3\
                   etot$!!ext 4 etotdel$!!ext 5 etotblost$!!ext 6 etotfladd$!!ext 7 etotinj$!!ext 8 etotrad$!!ext 9 etotdiff$!!ext  10\
		   m$!!ext   11 mdel$!!ext   12 mblost$!!ext   13 mfladd$!!ext   14 minj$!!ext   15 mrad$!!ext    16 mdiff$!!ext    17\
                   ie$!!ext  18 iedel$!!ext  19 ieblost$!!ext  20 iefladd$!!ext  21 ieinj$!!ext  22 ierad$!!ext   23 iediff$!!ext   24\
                   pe$!!ext 25 pedel$!!ext 26 peblost$!!ext 27 pefladd$!!ext 28 peinj$!!ext 29 perad$!!ext  30 pediff$!!ext  31\
		   ke$!!ext  32 kedel$!!ext  33 keblost$!!ext  34 kefladd$!!ext  35 keinj$!!ext  36 kerad$!!ext   37 kediff$!!ext   38\
                   mx1$!!ext 39 mx1del$!!ext 40 mx1blost$!!ext 41 mx1fladd$!!ext 42 mx1inj$!!ext 43 mx1rad$!!ext  44 mx1diff$!!ext  45\
                   mx2$!!ext 46 mx2del$!!ext 47 mx2blost$!!ext 48 mx2fladd$!!ext 49 mx2inj$!!ext 50 mx2rad$!!ext  51 mx2diff$!!ext  52\
		   mx3$!!ext 53 mx3del$!!ext 54 mx3blost$!!ext 55 mx3fladd$!!ext 56 mx3inj$!!ext 57 mx3rad$!!ext  58 mx3diff$!!ext  59\
		   ve$!!ext  60 vedel$!!ext  61 veblost$!!ext  62 vefladd$!!ext  63 veinj$!!ext  64 verad$!!ext   65 vediff$!!ext   66}
                   #
		  }
		  define enerreadtest (etotdel$!!ext[0])
                 }		 
                 if($fileversion==4){\
                  set minj$!!ext=0,0,1
                  set minj$!!ext[0]=-100
                  read {t$!!ext 1 cmode_amp$!!ext 2 smode_amp$!!ext 3\
                   etot$!!ext 4 etotdel$!!ext 5 etotblost$!!ext 6 etotfladd$!!ext 7 etotinj$!!ext 8 etotdiff$!!ext 9\
		   m$!!ext   10 mdel$!!ext   11 mblost$!!ext   12 mfladd$!!ext   13 minj$!!ext   14 mdiff$!!ext   15\
                   ie$!!ext  16 iedel$!!ext  17 ieblost$!!ext  18 iefladd$!!ext  19 ieinj$!!ext  20 iediff$!!ext  21\
                   pe$!!ext 22 pedel$!!ext 23 peblost$!!ext 24 pefladd$!!ext 25 peinj$!!ext 26 pediff$!!ext 27\
		   ke$!!ext  28 kedel$!!ext  29 keblost$!!ext  30 kefladd$!!ext  31 keinj$!!ext  32 kediff$!!ext  33\
                   mx1$!!ext 34 mx1del$!!ext 35 mx1blost$!!ext 36 mx1fladd$!!ext 37 mx1inj$!!ext 38 mx1diff$!!ext 39\
                   mx2$!!ext 40 mx2del$!!ext 41 mx2blost$!!ext 42 mx2fladd$!!ext 43 mx2inj$!!ext 44 mx2diff$!!ext 45\
		   mx3$!!ext 46 mx3del$!!ext 47 mx3blost$!!ext 48 mx3fladd$!!ext 49 mx3inj$!!ext 50 mx3diff$!!ext 51\
		   ve$!!ext  52 vedel$!!ext  53 veblost$!!ext  54 vefladd$!!ext  55 veinj$!!ext  56 vediff$!!ext  57}
                  define enerreadtest (minj$!!ext[0])
                 }
                 if($fileversion==3){\
                   set mfladd$!!ext=0,0,1
                   set mfladd$!!ext[0]=-100
                   read {t$!!ext 1 cmode_amp$!!ext 2 smode_amp$!!ext 3\
                   etot$!!ext 4 etotdel$!!ext 5 etotblost$!!ext 6 etotfladd$!!ext 7 etotdiff$!!ext 8\
		   m$!!ext    9 mdel$!!ext   10 mblost$!!ext   11 mfladd$!!ext   12 mdiff$!!ext   13\
                   ie$!!ext  14 iedel$!!ext  15 ieblost$!!ext  16 iefladd$!!ext  17 iediff$!!ext  18\
                   pe$!!ext 19 pedel$!!ext 20 peblost$!!ext 21 pefladd$!!ext 22 pediff$!!ext 23\
		   ke$!!ext  24 kedel$!!ext  25 keblost$!!ext  26 kefladd$!!ext  27 kediff$!!ext  28\
                   mx1$!!ext 29 mx1del$!!ext 30 mx1blost$!!ext 31 mx1fladd$!!ext 32 mx1diff$!!ext 33\
                   mx2$!!ext 34 mx2del$!!ext 35 mx2blost$!!ext 36 mx2fladd$!!ext 37 mx2diff$!!ext 38\
		   mx3$!!ext 39 mx3del$!!ext 40 mx3blost$!!ext 41 mx3fladd$!!ext 42 mx3diff$!!ext 43\
		   ve$!!ext  44 vedel$!!ext  45 veblost$!!ext  46 vefladd $!!ext 47 vediff$!!ext  48}
                   define enerreadtest (mfladd$!!ext[0])
                   #
                 }
                 # ZEUS HISTORY FILE (default)
                 if($fileversion==100){\
                   read {t 1 dt 2 m 3 ie 4 ke1 5 ke2 6 ke3 7 be1 8 be2 9 be3 10 re 11}
 		   set ke = (ke1+ke2+ke3)
		   set be = (be1+be2+be3)
 		   set etot = (ie+ke1+ke2+ke3+be1+be2+be3)
		   set etotdel = (etot-etot[0])
                   define enerreadtest (1)
                   #
                 }
                }
rdreads3   1    #
                if($filetype==$LOSSTYPE){\
                 if($fileversion==5){\
		  if($fullvec){\
                   read {t$!!ext 1 totm$!!ext 2 min1$!!ext 3 mout1$!!ext 4 min2$!!ext 5 mout2$!!ext 6\
			   toth$!!ext 7 hin1$!!ext 8 hout1$!!ext 9 hin2$!!ext 10 hout2$!!ext 11\
			   totpe$!!ext 12 pein1$!!ext 13 peout1$!!ext 14 pein2$!!ext 15 peout2$!!ext 16\
			   totke$!!ext 17 kein1$!!ext 18 keout1$!!ext 19 kein2$!!ext 20 keout2$!!ext 21\
			   totmx1$!!ext 22 mx1in1$!!ext 23 mx1out1$!!ext 24 mx1in2$!!ext 25 mx1out2$!!ext 26\
			   totmx2$!!ext 27 mx2in1$!!ext 28 mx2out1$!!ext 29 mx2in2$!!ext 30 mx2out2$!!ext 31\
			   totmx3$!!ext 32 mx3in1$!!ext 33 mx3out1$!!ext 34 mx3in2$!!ext 35 mx3out2$!!ext 36\
			   totbe$!!ext 37 bein1$!!ext 38 beout1$!!ext 39 bein2$!!ext 40 beout2$!!ext 41\
			   totbx1$!!ext 42 bx1in1$!!ext 43 bx1out1$!!ext 44 bx1in2$!!ext 45 bx1out2$!!ext 46\
			   totbx2$!!ext 47 bx2in1$!!ext 48 bx2out1$!!ext 49 bx2in2$!!ext 50 bx2out2$!!ext 51\
			   totbx3$!!ext 52 bx3in1$!!ext 53 bx3out1$!!ext 54 bx3in2$!!ext 55 bx3out2$!!ext 56\
			   totve$!!ext 57 vein1$!!ext 58 veout1$!!ext 59 vein2$!!ext 60 veout2$!!ext 61}
                   set tote$!!ext=toth$!!ext+totpe$!!ext+totke$!!ext+totbe$!!ext+totve$!!ext
                   set ein1$!!ext=hin1$!!ext+pein1$!!ext+kein1$!!ext+bein1$!!ext+vein1$!!ext
                   set eout1$!!ext=hout1$!!ext+peout1$!!ext+keout1$!!ext+beout1$!!ext+veout1$!!ext
                   set ein2$!!ext=hin2$!!ext+pein2$!!ext+kein2$!!ext+bein2$!!ext+vein2$!!ext
                   set eout2$!!ext=hout2$!!ext+peout2$!!ext+keout2$!!ext+beout2$!!ext+veout2$!!ext
		  }\
		  else{\
                   read {t$!!ext 1 totm$!!ext 2 min1$!!ext 3 mout1$!!ext 4 min2$!!ext 5 mout2$!!ext 6\
			   toth$!!ext 7 hin1$!!ext 8 hout1$!!ext 9 hin2$!!ext 10 hout2$!!ext 11\
			   totpe$!!ext 12 pein1$!!ext 13 peout1$!!ext 14 pein2$!!ext 15 peout2$!!ext 16\
			   totke$!!ext 17 kein1$!!ext 18 keout1$!!ext 19 kein2$!!ext 20 keout2$!!ext 21\
			   totmx1$!!ext 22 mx1in1$!!ext 23 mx1out1$!!ext 24 mx1in2$!!ext 25 mx1out2$!!ext 26\
			   totmx2$!!ext 27 mx2in1$!!ext 28 mx2out1$!!ext 29 mx2in2$!!ext 30 mx2out2$!!ext 31\
			   totmx3$!!ext 32 mx3in1$!!ext 33 mx3out1$!!ext 34 mx3in2$!!ext 35 mx3out2$!!ext 36\
			   totve$!!ext 37 vein1$!!ext 38 veout1$!!ext 39 vein2$!!ext 40 veout2$!!ext 41}
                   set tote$!!ext=toth$!!ext+totpe$!!ext+totke$!!ext+totve$!!ext
                   set ein1$!!ext=hin1$!!ext+pein1$!!ext+kein1$!!ext+vein1$!!ext
                   set eout1$!!ext=hout1$!!ext+peout1$!!ext+keout1$!!ext+veout1$!!ext
                   set ein2$!!ext=hin2$!!ext+pein2$!!ext+kein2$!!ext+vein2$!!ext
                   set eout2$!!ext=hout2$!!ext+peout2$!!ext+keout2$!!ext+veout2$!!ext
		  }
		}
		if($fileversion==6){\
		 if($LOOPTYPE==1){\
                  read {t$!!ext 1 totm$!!ext 2 min1$!!ext 3 mout1$!!ext 4 min2$!!ext 5 mout2$!!ext 6 min3$!!ext 7 mout3$!!ext 8\
		       toth$!!ext 9 hin1$!!ext 10 hout1$!!ext 11 hin2$!!ext 12 hout2$!!ext 13  hin3$!!ext 14 hout3$!!ext 15\
		       totpe$!!ext 16 pein1$!!ext 17 peout1$!!ext 18 pein2$!!ext 19 peout2$!!ext 20 pein3$!!ext 21 peout3$!!ext 22\
		       totke$!!ext 23 kein1$!!ext 24 keout1$!!ext 25 kein2$!!ext 26 keout2$!!ext 27 kein3$!!ext 28 keout3$!!ext 29\
		       totmx1$!!ext 30 mx1in1$!!ext 31 mx1out1$!!ext 32 mx1in2$!!ext 33 mx1out2$!!ext 34 mx1in3$!!ext 35 mx1out3$!!ext 36\
		       totmx2$!!ext 37 mx2in1$!!ext 38 mx2out1$!!ext 39 mx2in2$!!ext 40 mx2out2$!!ext 41 mx2in3$!!ext 42 mx2out3$!!ext 43\
		       totmx3$!!ext 44 mx3in1$!!ext 45 mx3out1$!!ext 46 mx3in2$!!ext 47 mx3out2$!!ext 48 mx3in3$!!ext 49 mx3out3$!!ext 50\
		       totbe$!!ext 51 bein1$!!ext 52 beout1$!!ext 53 bein2$!!ext 54 beout2$!!ext 55 bein3$!!ext 56 beout3$!!ext 57\
		       totbx1$!!ext 58 bx1in1$!!ext 59 bx1out1$!!ext 60 bx1in2$!!ext 61 bx1out2$!!ext 62 bx1in3$!!ext 63 bx1out3$!!ext 64\
		       totbx2$!!ext 65 bx2in1$!!ext 66 bx2out1$!!ext 67 bx2in2$!!ext 68 bx2out2$!!ext 69 bx2in3$!!ext 70 bx2out3$!!ext 71\
		       totbx3$!!ext 72 bx3in1$!!ext 73 bx3out1$!!ext 74 bx3in2$!!ext 75 bx3out2$!!ext 76 bx3in3$!!ext 77 bx3out3$!!ext 78\
		       totve$!!ext 79 vein1$!!ext 80 veout1$!!ext 81 vein2$!!ext 82 veout2$!!ext 83 vein3$!!ext 84 veout3$!!ext 85}
                        set tote$!!ext=toth$!!ext+totpe$!!ext+totke$!!ext+totbe$!!ext+totve$!!ext
                        set ein1$!!ext=hin1$!!ext+pein1$!!ext+kein1$!!ext+bein1$!!ext+vein1$!!ext
                        set eout1$!!ext=hout1$!!ext+peout1$!!ext+keout1$!!ext+beout1$!!ext+veout1$!!ext
                        set ein2$!!ext=hin2$!!ext+pein2$!!ext+kein2$!!ext+bein2$!!ext+vein2$!!ext
                        set eout2$!!ext=hout2$!!ext+peout2$!!ext+keout2$!!ext+beout2$!!ext+veout2$!!ext
		        set ein3$!!ext=hin3$!!ext+pein3$!!ext+kein3$!!ext+bein3$!!ext+vein3$!!ext
                        set eout3$!!ext=hout3$!!ext+peout3$!!ext+keout3$!!ext+beout3$!!ext+veout3$!!ext
		  }
		  if($LOOPTYPE>1){\
                   read {t$!!ext 1 min$!!ext 2 mout$!!ext 3\
		                   hin$!!ext 4 hout$!!ext 5\
		                   pein$!!ext 6 peout$!!ext 7\
		                   kein$!!ext 8 keout$!!ext 9\
		                   mx1in$!!ext 10 mx1out$!!ext 11\
		                   mx2in$!!ext 12 mx2out$!!ext 13\
		                   mx3in$!!ext 14 mx3out$!!ext 15\
		                   bein$!!ext 16 beout$!!ext 17\
		                   bx1in$!!ext 18 bx1out$!!ext 19\
		                   bx2in$!!ext 20 bx2out$!!ext 21\
		                   bx3in$!!ext 22 bx3out$!!ext 23\
		                   vein$!!ext 24 veout$!!ext 25}
		   set mflux$!!ext = min$!!ext + mout$!!ext
		   set hflux$!!ext = hin$!!ext + hout$!!ext
		   set peflux$!!ext = pein$!!ext + peout$!!ext
		   set keflux$!!ext = kein$!!ext + keout$!!ext
		   set mx1flux$!!ext = mx1in$!!ext + mx1out$!!ext
		   set mx2flux$!!ext = mx2in$!!ext + mx2out$!!ext
		   set mx3flux$!!ext = mx3in$!!ext + mx3out$!!ext
		   set beflux$!!ext = bein$!!ext + beout$!!ext
		   set bx1flux$!!ext = bx1in$!!ext + bx1out$!!ext
		   set bx2flux$!!ext = bx2in$!!ext + bx2out$!!ext
		   set bx3flux$!!ext = bx3in$!!ext + bx3out$!!ext
		   set veflux$!!ext = vein$!!ext + veout$!!ext
	 	
	 	   set ein$!!ext=hin$!!ext+pein$!!ext+kein$!!ext+bein$!!ext+vein$!!ext
		   set eout$!!ext=hout$!!ext+peout$!!ext+keout$!!ext+beout$!!ext+veout$!!ext
	 	   set eflux$!!ext=ein$!!ext+eout$!!ext
		 }
		}
	       }
rdreads31  1   #
	       if($filetype==$LOSSTYPE){\
		if($fileversion==7){\
		 if($LOOPTYPE==1){\
                  read {t$!!ext 1 totm$!!ext 2 min1$!!ext 3 mout1$!!ext 4 min2$!!ext 5 mout2$!!ext 6 min3$!!ext 7 mout3$!!ext 8\
		       toth$!!ext 9 hin1$!!ext 10 hout1$!!ext 11 hin2$!!ext 12 hout2$!!ext 13  hin3$!!ext 14 hout3$!!ext 15\
		       totpe$!!ext 16 pein1$!!ext 17 peout1$!!ext 18 pein2$!!ext 19 peout2$!!ext 20 pein3$!!ext 21 peout3$!!ext 22\
		       totke$!!ext 23 kein1$!!ext 24 keout1$!!ext 25 kein2$!!ext 26 keout2$!!ext 27 kein3$!!ext 28 keout3$!!ext 29\
		       totmx1$!!ext 30 mx1in1$!!ext 31 mx1out1$!!ext 32 mx1in2$!!ext 33 mx1out2$!!ext 34 mx1in3$!!ext 35 mx1out3$!!ext 36\
		       totmx2$!!ext 37 mx2in1$!!ext 38 mx2out1$!!ext 39 mx2in2$!!ext 40 mx2out2$!!ext 41 mx2in3$!!ext 42 mx2out3$!!ext 43\
		       totmx3$!!ext 44 mx3in1$!!ext 45 mx3out1$!!ext 46 mx3in2$!!ext 47 mx3out2$!!ext 48 mx3in3$!!ext 49 mx3out3$!!ext 50\
		       totbe$!!ext 51 bein1$!!ext 52 beout1$!!ext 53 bein2$!!ext 54 beout2$!!ext 55 bein3$!!ext 56 beout3$!!ext 57\
		       totbx1$!!ext 58 bx1in1$!!ext 59 bx1out1$!!ext 60 bx1in2$!!ext 61 bx1out2$!!ext 62 bx1in3$!!ext 63 bx1out3$!!ext 64\
		       totbx2$!!ext 65 bx2in1$!!ext 66 bx2out1$!!ext 67 bx2in2$!!ext 68 bx2out2$!!ext 69 bx2in3$!!ext 70 bx2out3$!!ext 71\
		       totbx3$!!ext 72 bx3in1$!!ext 73 bx3out1$!!ext 74 bx3in2$!!ext 75 bx3out2$!!ext 76 bx3in3$!!ext 77 bx3out3$!!ext 78\
		       totve$!!ext 79 vein1$!!ext 80 veout1$!!ext 81 vein2$!!ext 82 veout2$!!ext 83 vein3$!!ext 84 veout3$!!ext 85\
		       totamx1$!!ext 86 amx1in1$!!ext 87 amx1out1$!!ext 88 amx1in2$!!ext 89 amx1out2$!!ext 90 amx1in3$!!ext 91 amx1out3$!!ext 92\
		       totamx2$!!ext 93 amx2in1$!!ext 94 amx2out1$!!ext 95 amx2in2$!!ext 96 amx2out2$!!ext 97 amx2in3$!!ext 98 amx2out3$!!ext 99\
		       totamx3$!!ext 100 amx3in1$!!ext 101 amx3out1$!!ext 102 amx3in2$!!ext 103 amx3out2$!!ext 104 amx3in3$!!ext 105 amx3out3$!!ext 106}
                        set tote$!!ext=toth$!!ext+totpe$!!ext+totke$!!ext+totbe$!!ext+totve$!!ext
                        set ein1$!!ext=hin1$!!ext+pein1$!!ext+kein1$!!ext+bein1$!!ext+vein1$!!ext
                        set eout1$!!ext=hout1$!!ext+peout1$!!ext+keout1$!!ext+beout1$!!ext+veout1$!!ext
                        set ein2$!!ext=hin2$!!ext+pein2$!!ext+kein2$!!ext+bein2$!!ext+vein2$!!ext
                        set eout2$!!ext=hout2$!!ext+peout2$!!ext+keout2$!!ext+beout2$!!ext+veout2$!!ext
		        set ein3$!!ext=hin3$!!ext+pein3$!!ext+kein3$!!ext+bein3$!!ext+vein3$!!ext
                        set eout3$!!ext=hout3$!!ext+peout3$!!ext+keout3$!!ext+beout3$!!ext+veout3$!!ext
		  }
		  if($LOOPTYPE>1){\
                   read {t$!!ext 1 min$!!ext 2 mout$!!ext 3\
		                   hin$!!ext 4 hout$!!ext 5\
		                   pein$!!ext 6 peout$!!ext 7\
		                   kein$!!ext 8 keout$!!ext 9\
		                   mx1in$!!ext 10 mx1out$!!ext 11\
		                   mx2in$!!ext 12 mx2out$!!ext 13\
		                   mx3in$!!ext 14 mx3out$!!ext 15\
		                   bein$!!ext 16 beout$!!ext 17\
		                   bx1in$!!ext 18 bx1out$!!ext 19\
		                   bx2in$!!ext 20 bx2out$!!ext 21\
		                   bx3in$!!ext 22 bx3out$!!ext 23\
		                   vein$!!ext 24 veout$!!ext 25\
		                   amx1in$!!ext 26 amx1out$!!ext 27\
		                   amx2in$!!ext 28 amx2out$!!ext 29\
		                   amx3in$!!ext 30 amx3out$!!ext 31}
		   set mflux$!!ext = min$!!ext + mout$!!ext
		   set hflux$!!ext = hin$!!ext + hout$!!ext
		   set peflux$!!ext = pein$!!ext + peout$!!ext
		   set keflux$!!ext = kein$!!ext + keout$!!ext
		   set mx1flux$!!ext = mx1in$!!ext + mx1out$!!ext
		   set mx2flux$!!ext = mx2in$!!ext + mx2out$!!ext
		   set mx3flux$!!ext = mx3in$!!ext + mx3out$!!ext
		   set beflux$!!ext = bein$!!ext + beout$!!ext
		   set bx1flux$!!ext = bx1in$!!ext + bx1out$!!ext
		   set bx2flux$!!ext = bx2in$!!ext + bx2out$!!ext
		   set bx3flux$!!ext = bx3in$!!ext + bx3out$!!ext
		   set veflux$!!ext = vein$!!ext + veout$!!ext
		   set amx1flux$!!ext = amx1in$!!ext + amx1out$!!ext
		   set amx2flux$!!ext = amx2in$!!ext + amx2out$!!ext
		   set amx3flux$!!ext = amx3in$!!ext + amx3out$!!ext
	 	
	 	   set ein$!!ext=hin$!!ext+pein$!!ext+kein$!!ext+bein$!!ext+vein$!!ext
		   set eout$!!ext=hout$!!ext+peout$!!ext+keout$!!ext+beout$!!ext+veout$!!ext
	 	   set eflux$!!ext=ein$!!ext+eout$!!ext
		}
	       }
              }
rdreads32   1 # second part	
                if($filetype==$SPTYPE){\
                 if($fileversion==1){\
                  read {t$!!ext 1 x1l$!!ext 2 x1h$!!ext 3  x2l$!!ext 4 x2h$!!ext 5  totl$!!ext 6 toth$!!ext 7}
                 }
                }
                if($filetype==$TSTYPE){\
                 if($fileversion==1){\
                  read {t 1 \
			p0c2$!!ext  2 p0c3$!!ext  3 p0c4$!!ext  4 p0c5$!!ext  5 p0c6$!!ext  6 p0c7$!!ext  7 p0c8$!!ext  8 p0c9$!!ext  9 p0c10$!!ext 10 \
		        p1c2$!!ext 11 p1c3$!!ext 12 p1c4$!!ext 13 p1c5$!!ext 14 p1c6$!!ext 15 p1c7$!!ext 16 p1c8$!!ext 17 p1c9$!!ext 18 p1c10$!!ext 19 \
		        p2c2$!!ext 20 p2c3$!!ext 21 p2c4$!!ext 22 p2c5$!!ext 23 p2c6$!!ext 24 p2c7$!!ext 25 p2c8$!!ext 26 p2c9$!!ext 27 p2c10$!!ext 28 \
		        p3c2$!!ext 29 p3c3$!!ext 30 p3c4$!!ext 31 p3c5$!!ext 32 p3c6$!!ext 33 p3c7$!!ext 34 p3c8$!!ext 35 p3c9$!!ext 36 p3c10$!!ext 37 \
		        p4c2$!!ext 38 p4c3$!!ext 39 p4c4$!!ext 40 p4c5$!!ext 41 p4c6$!!ext 42 p4c7$!!ext 43 p4c8$!!ext 44 p4c9$!!ext 45 p4c10$!!ext 46}
                 }
                }
                if($filetype==$LOGDTTYPE){\
                 if($fileversion==1){\
                  read {t$!!ext 1 dt$!!ext 2 r$!!ext 3 w$!!ext 4 l$!!ext 5 k$!!ext 6 j$!!ext 7 i$!!ext 8 x1$!!ext 9 x2$!!ext 10 \
                   cs_dt$!!ext 11 cs_v$!!ext 12 x1v_dt$!!ext 13 x1v_v$!!ext 14 x2v_dt$!!ext 15 x2v_v$!!ext 16 bv_dt$!!ext 17 bv_v$!!ext 18 lv_dt$!!ext 19\
                   lv_dv$!!ext 20 vx1_dt$!!ext 21 vx1_dv$!!ext 22 vx2_dt$!!ext 23 vx2_dv$!!ext 24 rvx1_dt$!!ext 25 rvx1_nu$!!ext 26 rvx2_dt$!!ext 27\
                   rvx2_nu$!!ext 28 res_dt$!!ext 29 res_v$!!ext 30}
                 }
                }
                if($filetype==$STEPTYPE){\
                 if($fileversion==1){\
                  read {n$!!ext 1 nn$!!ext 2 nsup$!!ext 3 nsub$!!ext 4 t$!!ext 5 dt$!!ext 6 upto$!!ext 7 i$!!ext 8 N$!!ext 9}
                 }
                }
                if($filetype==$PERFTYPE){\
                 if($fileversion==3){\
                  read {t$!!ext 1 ete$!!ext 2 n$!!ext 3 wt$!!ext 4 zc$!!ext 5 tuphr$!!ext 6 lete$!!ext 7 ln$!!ext 8 lwt$!!ext 9 lzc$!!ext 10 ltuphr$!!ext 11}
                  #t, ete, n, wt, zc, tu/hr,  lete, ln, lwt, lzc, ltu/hr 
                 }
                 if($fileversion==2){\
                  read {proc$!!ext 1 t$!!ext 2 ete$!!ext 3 n$!!ext 4 wt$!!ext 5 zc$!!ext 6 tuphr$!!ext 7 lete$!!ext 8 ln$!!ext 9 lwt$!!ext 10 lzc$!!ext 11 ltuphr$!!ext 12}
                 }
                 if($fileversion==1){\
                  read {proc$!!ext 1 ete$!!ext 2 n$!!ext 3 wt$!!ext 4 zc$!!ext 5 lete$!!ext 6 ln$!!ext 7 lwt$!!ext 8 lzc$!!ext 9}
                 }
                }
                #
rdreads4 1      #
		if($filetype==$MODETYPE){\
                 # only version so far!
		 # 10 modes, 0-9
		 if($fileversion==1){\
                        read {t$!!ext 1 m0$!!ext 2 m1$!!ext 3 m2$!!ext 4 m3$!!ext 5 m4$!!ext 6 m5$!!ext 7 m6$!!ext 8 m7$!!ext 9 m8$!!ext 10 m9$!!ext 11}
                 }
                 if($fileversion==2){\
		        read {t$!!ext 1 \
		        rm0$!!ext  2 rm1$!!ext 3 rm2$!!ext 4 rm3$!!ext 5 rm4$!!ext 6 rm5$!!ext 7 rm6$!!ext 8 rm7$!!ext 9 rm8$!!ext 10 rm9$!!ext 11 \
		        em0$!!ext  12 em1$!!ext 13 em2$!!ext 14 em3$!!ext 15 em4$!!ext 16 em5$!!ext 17 em6$!!ext 18 em7$!!ext 19 em8$!!ext 20 em9$!!ext 21 \
		        pem0$!!ext 22 pem1$!!ext 23 pem2$!!ext 24 pem3$!!ext 25 pem4$!!ext 26 pem5$!!ext 27 pem6$!!ext 28 pem7$!!ext 29 pem8$!!ext 30 pem9$!!ext 31 \
		        vm0$!!extx 32 vm1$!!extx 33 vm2$!!extx 34 vm3$!!extx 35 vm4$!!extx 36 vm5$!!extx 37 vm6$!!extx 38 vm7$!!extx 39 vm8$!!extx 40 vm9$!!extx 41 \
		        vm0$!!exty 42 vm1$!!exty 43 vm2$!!exty 44 vm3$!!exty 45 vm4$!!exty 46 vm5$!!exty 47 vm6$!!exty 48 vm7$!!exty 49 vm8$!!exty 50 vm9$!!exty 51 \
		        vm0$!!extz 52 vm1$!!extz 53 vm2$!!extz 54 vm3$!!extz 55 vm4$!!extz 56 vm5$!!extz 57 vm6$!!extz 58 vm7$!!extz 59 vm8$!!extz 60 vm9$!!extz 61 \
		        bm0$!!extx 62 bm1$!!extx 63 bm2$!!extx 64 bm3$!!extx 65 bm4$!!extx 66 bm5$!!extx 67 bm6$!!extx 68 bm7$!!extx 69 bm8$!!extx 70 bm9$!!extx 71 \
		        bm0$!!exty 72 bm1$!!exty 73 bm2$!!exty 74 bm3$!!exty 75 bm4$!!exty 76 bm5$!!exty 77 bm6$!!exty 78 bm7$!!exty 79 bm8$!!exty 80 bm9$!!exty 81 \
		        bm0$!!extz 82 bm1$!!extz 83 bm2$!!extz 84 bm3$!!extz 85 bm4$!!extz 86 bm5$!!extz 87 bm6$!!extz 88 bm7$!!extz 89 bm8$!!extz 90 bm9$!!extz 91}
		 }
		}
		#
rdcalcs1 1      #
                # any cross-grid stuff below assumes DUMPSM==1 in code
                #
                # if $gam changes, $gam stuff is wrong!
                if(($filetype==$DTYPE)||($filetype==$PDTYPE)||($filetype==$ADTYPE)||($filetype==$FLTYPE)||($filetype==$AVG2DTYPE)||($filetype==$AVG1DTYPE)){\
                 set p$!!ext = ($gam-1.)*en$!!ext
                 if(wgam){\
                   set csq$!!ext = $gam*($gam-1.)*en$!!ext/r$!!ext
	           set cs$!!ext = SQRT(csq$!!ext)
                   set cstimex1$!!ext = cs$!!ext/dx1
                   set cstimex2$!!ext = cs$!!ext/dx2
                 }
 		 if(_mag){\
                  set b2$!!ext = bx**2+by**2+bz**2
		  set lgb2$!!ext=LG(b2$!!ext)
                  set absb$!!ext = sqrt(b2$!!ext)
                  set va$!!ext = absb/sqrt(r$!!ext)
		  # if using alfven limiter
		  set vareal$!!ext = absb/sqrt(r$!!ext+b2$!!ext*invsol2)
		  # if CS limiter
		  set csreal$!!ext = $gam*($gam-1)*en$!!ext/(r$!!ext+$gam*($gam-1)*en$!!ext*invsol2)
		  #
		  set va$!!extx = ABS(bx)/sqrt(r$!!ext)
		  set vatimex1$!!ext = va$!!ext/dx1
		  set vatimex2$!!ext = va$!!ext/dx2
		  set magbeta$!!ext = ($gam-1)*(en$!!ext)/((b2$!!ext*0.5+1E-9))
		  set imagbeta$!!ext =((b2$!!ext*0.5))/( ($gam-1)*(en$!!ext))
		 }\
		 else{\
			set b2$!!ext = 0
			set absb$!!ext = 0
			set va$!!ext = 0
		        set va$!!extx = 0
			set vatimex1$!!ext = 0
			set vatimex2$!!ext = 0
		 }
                 set tottimex1$!!ext = (ABS(v$!!extx)+va$!!ext+cs$!!ext)/dx1
		 set tottimex2$!!ext = (ABS(v$!!exty)+va$!!ext+cs$!!ext)/(x12*dx2)
                 set ltottimex1$!!ext = (ABS(v$!!extx)+sqrt(1/(1+1/va$!!ext)+1/(1+1/cs$!!ext)))/dx1
		 set ltottimex2$!!ext = (ABS(v$!!exty)+sqrt(1/(1+1/va$!!ext)+1/(1+1/cs$!!ext)))/(x12*dx2)
                 set vtotsq$!!ext = v$!!extx*v$!!extx+v$!!exty*v$!!exty+v$!!extz*v$!!extz
                 set vtot$!!ext = SQRT(vtotsq$!!ext)
                 set enthalpy$!!ext = $gam*en$!!ext/r$!!ext
		 if(ABS($gam-1.)<.000001){\
					  set entropy$!!ext = r$!!ext/r$!!ext
		                         }\
   		                         else{\
		                          # entropy per degree of freedom
		                          set entropy$!!ext = ($gam-1.)*en$!!ext/r$!!ext**($gam)
 		                          set me=9.1*10**(-28)
		                          set mp=1.67265*10**(-24)
		                          set hbar=1.05459*10**(-27)
		                          set kb=1.38066*10**(-16)
		                          set eentropy=entropy$!!ext*(me**($gam+1))/(2*pi*kb*hbar**2)
		                          set pentropy=entropy$!!ext*(mp**($gam+1))/(2*pi*kb*hbar**2)
		                         }
                 set Mach$!!extx = (v$!!extx/(cs$!!ext+va$!!ext))
		 set Machcs$!!extx = (v$!!extx/(cs$!!ext))
		 set Machva$!!extx = (v$!!extx/(va$!!ext))
		 set Machvareal$!!extx = (v$!!extx/(vareal$!!ext))
		 set Machreal$!!extx = (v$!!extx/(vareal$!!ext+csreal$!!ext))
                 set Mach$!!exty = (v$!!exty/(cs$!!ext+va$!!ext))
                 set Mach$!!extz = (v$!!extz/(cs$!!ext+va$!!ext))
                 set Machv$!!ext = (vtot$!!ext/cs$!!ext)
                 set dl$!!ext = r$!!ext*(v$!!extz*x12*g42-SQRT(x12*g42))
		 if(_mag==1){\
		        set Be$!!ext = vtotsq$!!ext*0.5+enthalpy$!!ext+pot$!!ext+0.5*b2$!!ext
		     }\
		            else{\
                        set Be$!!ext = vtotsq$!!ext*0.5+enthalpy$!!ext+pot$!!ext
		     }
                 set se$!!ext = en$!!ext/r$!!ext
		 set ek$!!ext = 0.5*r$!!ext*vtotsq$!!ext
		 set p$!!extx = v$!!extx*r$!!ext
		 set p$!!exty = v$!!exty*r$!!ext
		 set lr$!!ext = LG(r$!!ext)
                 set len$!!ext = LG(en$!!ext)
                 set rv$!!extx = r$!!ext*v$!!extx
                 set rvxsq$!!ext = r$!!ext*v$!!extx*v$!!extx
		 set rv$!!exty = r$!!ext*v$!!exty
		 set rvysq$!!ext = r$!!ext*v$!!exty*v$!!exty
		 set rv$!!extz = r$!!ext*v$!!extz
		 set rvzsq$!!ext = r$!!ext*v$!!extz*v$!!extz
		 set mdot=4*pi*x12*r$!!ext*v$!!extx
		 # temporary
		 #set it=ABS(mdot)
		 #thetaphiavg PI/6 it mdota newx1
		 if($coord==3){\
		        set angmom3$!!ext = x12*g42*r$!!ext*v$!!extz
		     }
		 if($coord==1){\
		            set radius=sqrt(x1**2+x2**2+x3**2)
		            set angmom3=radius*sqrt(vx**2+vy**2)*r*dx1*dx2*dx3
		            set omega3cart=sqrt(vx**2+vy**2)/radius
		         }
                 # spherical rotation
		 set omega3$!!ext = v$!!extz/(x12*g42)
                 set omegak$!!ext = sqrt(1.0/(x12**3 * (1.0-2/x12)**2 ) )
		 if(_mag){\
                  set dx1mri=3*dx1*(ABS(omega3))/(2*PI*ABS(va)/sqrt(3))
                  set dx2mri=3*x1*(ABS(omega3))*dx2/(2*PI*ABS(va)/sqrt(3))
		 }
                 if($npdone){\
                  # visc times for various perscriptions
                  #need to divide by alpha
                  # for IGU
                  #set visctime$!!ext=x12*x12*omegak$!!ext/(cs$!!ext*cs$!!ext)/$alphareal0
		  set lvisctimex$!!ext=dx12*dx12/(nuvisc$!!ext)
		  set lvisctimey$!!ext=g32*g32*dx22*dx22/(nuvisc$!!ext)
		  set visctime$!!ext=x12*x12/(nuvisc$!!ext)
                 }
                 set mass=r$!!ext*dvl1*dvl2

                 # Mass Flux density
                 set Fmd$!!extx=r$!!ext*v$!!extx
                 set Fmd$!!exty=r$!!ext*v$!!exty
                 set Fmd$!!extz=r$!!ext*v$!!extz

                 # Mass Flux(full reduced derivative argument, all left is ds)
                 set Fm$!!extx=g22*g32*g42*r$!!ext*v$!!extx
                 set Fm$!!exty=g32*g42*r$!!ext*v$!!exty
                 set Fm$!!extz=g32*r$!!ext*v$!!extz

                 # Ang mom Reynolds flux density(argument of divergence)
                 set FPvzd$!!extx=g32*g42*(r$!!ext*v$!!extx*v$!!extz)
                 set FPvzd$!!exty=g32*g42*(r$!!ext*v$!!exty*v$!!extz)
                 set FPvzd$!!extz=g32*g42*(r$!!ext*v$!!extz*v$!!extz)

                 # Ang mom Reynolds flux(full reduced derivative argument)
                 set FPvz$!!extx=g22*g32*g32*g42*g42*(r$!!ext*v$!!extx*v$!!extz)
                 set FPvz$!!exty=g22*g32*g42*g42*(r$!!ext*v$!!exty*v$!!extz)
                 set FPvz$!!extz=g22*g32*g42*(r$!!ext*v$!!extz*v$!!extz)

                 # Energy Flux (density(argument) and full derivative argument)
                 set temp$!!ext=0.5*vtotsq$!!ext+enthalpy$!!ext+pot$!!ext

                 set FPEd$!!extx =  temp$!!ext*Fmd$!!extx
                 set FPEd$!!exty =  temp$!!ext*Fmd$!!exty
                 set FPEd$!!extz =  temp$!!ext*Fmd$!!extz

		 set FPE$!!extx =  temp$!!ext*Fm$!!extx
                 set FPE$!!exty =  temp$!!ext*Fm$!!exty
                 set FPE$!!extz =  temp$!!ext*Fm$!!extz

		 if(_mag==1){\
		  set temp$!!ext = vx*bx+vy*by+vz*bz
		}\
		else{\
		         set temp$!!ext = 0
		      }
		  set FME1d$!!extx = -temp$!!ext*bx
		  set FME1d$!!exty = -temp$!!ext*by
		  set FME1d$!!extz = -temp$!!ext*bz

		  set temp$!!ext = bx*bx+by*by+bz*bz
		  set FME2d$!!extx = temp$!!ext*vx
		  set FME2d$!!exty = temp$!!ext*vy
		  set FME2d$!!extz = temp$!!ext*vz

		  set FMEd$!!extx=FME1d$!!extx+FME2d$!!extx
		  set FMEd$!!exty=FME1d$!!exty+FME2d$!!exty
		  set FMEd$!!extz=FME1d$!!extz+FME2d$!!extz

		  set FME1$!!extx=g22*g32*g42*FME1d$!!extx
		  set FME1$!!exty=g32*g42*FME1d$!!exty
		  set FME1$!!extz=g32*FME1d$!!extz

		  set FME2$!!extx=g22*g32*g42*FME2d$!!extx
		  set FME2$!!exty=g32*g42*FME2d$!!exty
		  set FME2$!!extz=g32*FME2d$!!extz

		  set FME$!!extx=FME1$!!extx+FME2$!!extx
		  set FME$!!exty=FME1$!!exty+FME2$!!exty
		  set FME$!!extz=FME1$!!extz+FME2$!!extz
                }
rdcalcs2   1    #
                # viscosity and total only parts
                if(($filetype==$DTYPE)||($filetype==$PDTYPE)||($filetype==$ADTYPE)||($filetype==$FLTYPE)||($filetype==$AVG2DTYPE)||($filetype==$AVG1DTYPE)){\
		 if($npdone){\
                  set alphavisc=sig13$!!ext/p
                  # Ang mom viscous flux density(argument of divergence)
                  set FVvzd$!!extx=g32*g42*sig13$!!ext
                  set FVvzd$!!exty=g32*g42*sig23$!!ext
                  set FVvzd$!!extz=g32*g42*sig33$!!ext

                  # Ang mom viscous flux(full reduced derivative argument)
                  set FVvz$!!extx=g22*g32*g22*g42*g42*sig13$!!ext
                  set FVvz$!!exty=g22*g32*g42*g42*sig23$!!ext
                  set FVvz$!!extz=g22*g32*g42*sig33$!!ext
		 }\
                 else{\
                  set FVvzdx=0
                  set FVvzdy=0
                  set FVvzdz=0
                  set FVvzx=0
                  set FVvzy=0
                  set FVvzz=0
		 }
                 # Total ang mom flux density(argument of divergence)
                 set Fvzd$!!extx = FPvzd$!!extx+FVvzd$!!extx
                 set Fvzd$!!exty = FPvzd$!!exty+FVvzd$!!exty
                 set Fvzd$!!extz = FPvzd$!!extz+FVvzd$!!extz

                 # Total ang mom flux(full reduced derivative argument)
                 set Fvz$!!extx = FPvz$!!extx+FVvz$!!extx
                 set Fvz$!!exty = FPvz$!!exty+FVvz$!!exty
                 set Fvz$!!extz = FPvz$!!extz+FVvz$!!extz
                 #
                 #
		 if($npdone){\
                  # argument to divergence only
                  set FVEd$!!extx = sig13$!!ext*v$!!extz+sig12$!!ext*v$!!exty+sig11$!!ext*v$!!extx
                  set FVEd$!!exty = sig23$!!ext*v$!!extz+sig22$!!ext*v$!!exty+sig21$!!ext*v$!!extx
                  set FVEd$!!extz = sig33$!!ext*v$!!extz+sig32$!!ext*v$!!exty+sig31$!!ext*v$!!extx
		 }\
		 else{\
                  set FVEdx=0
                  set FVEdy=0
                  set FVEdz=0
		 }
		 if($npdone){\
                  # full reduced arguments
                  set FVE$!!extx = g22*g32*g42*FVEd$!!extx
                  set FVE$!!exty = g32*g42*FVEd$!!exty
		  set FVE$!!extz = g32*FVEd$!!extz
		 }\
		 else{\
                  set FVEx=0
                  set FVEy=0
                  set FVEz=0
		 }

		 # total energy flux
                 set FEd$!!extx = FVEd$!!extx+FPEd$!!extx+FMEd$!!extx
                 set FEd$!!exty = FVEd$!!exty+FPEd$!!exty+FMEd$!!exty
                 set FEd$!!extz = FVEd$!!extz+FPEd$!!extz+FMEd$!!extz
 
                 set FE$!!extx = FVE$!!extx+FPE$!!extx+FME$!!extx
                 set FE$!!exty = FVE$!!exty+FPE$!!exty+FME$!!exty
                 set FE$!!extz = FVE$!!extz+FPE$!!extz+FME$!!extz
                
		 if($npdone){\
                  # NP calcs (assume npdump already read in once dump hit)
                  if(($filetype==$DTYPE)||($filetype==$AVG2DTYPE)||($filetype==$AVG1DTYPE)){\
		   if($fileversion==1){\
                    set tempvisc=nuvisc if(nuvisc>0)
                    if(dimen(tempvisc)>0){\
                     set Phi$!!ext = (sig11$!!ext*sig11$!!ext+sig22$!!ext*sig22$!!ext+sig33$!!ext*sig33$!!ext+2.0*sig13$!!ext*sig13$!!ext+2.0*sig23$!!ext*sig23$!!ext+2.0*sig21$!!ext*sig21$!!ext)/(2.0*r$!!ext*nuvisc$!!ext)
                    }
                   }
		 
 		   if($fileversion==0){\
		    # probably MG prescription if this is the case, but may be turned off at a certain time, but sigmas will take care of that.
		    set nuvisc$!!ext=_alphareal0*cs$!!ext*cs$!!ext/omegak$!!ext*g42**(1.5)
		    set Phi$!!ext = (sig11$!!ext*sig11$!!ext+sig22$!!ext*sig22$!!ext+sig33$!!ext*sig33$!!ext+2.0*sig13$!!ext*sig13$!!ext+2.0*sig23$!!ext*sig23$!!ext+2.0*sig21$!!ext*sig21$!!ext)/(2.0*r$!!ext*nuvisc$!!ext)
                   }
                  }
		 }
		}
                # AVG2D calcs
                if($filetype==$AVG2DTYPE){\
                     #                done above essentially
                }
                # AVG1D calcs
                if($filetype==$AVG1DTYPE){\
                 set lx1$!!ext = LG(x1$!!ext)
                 #set lx2$!!ext = LG(x2$!!ext)
                 set lr1d1$!!ext = LG(r1d1$!!ext)
                 set len1d1$!!ext = LG(en1d1$!!ext)
                 set lbe1d1$!!ext = LG(ABS(be1d1$!!ext))
                 set lcsq1d1$!!ext = LG(csq1d1$!!ext)
                 set ls1d1$!!ext = LG(ABS(s1d1$!!ext))
                 set lvx1d1$!!ext = LG(v1d1$!!extx)
                 set lvy1d1$!!ext = LG(ABS(v1d1$!!exty))
                 set lvz1d1$!!ext = LG(v1d1$!!extz)
                 set lr1d2$!!ext = LG(r1d2$!!ext)
                 set len1d2$!!ext = LG(en1d2$!!ext)
                 set lbe1d2$!!ext = LG(ABS(be1d2$!!ext))
                 set lcsq1d2$!!ext = LG(csq1d2$!!ext)
                 set ls1d2$!!ext = LG(ABS(s1d2$!!ext))
                 set lvx1d2$!!ext = LG(v1d2$!!extx)
                 set lvy1d2$!!ext = LG(ABS(v1d2$!!exty))
                 set lvz1d2$!!ext = LG(v1d2$!!extz)
                 set Mach1$!!extx  = v1d1$!!extx/SQRT(csq1d1$!!ext)
                 set Mach2$!!extx  = v1d2$!!extx/SQRT(csq1d2$!!ext)
                }
                # CALC calcs
                if($filetype==$CALCTYPE){\
                 set L$!!extx  = LR$!!extx+LV$!!extx
                 set L$!!exty  = LR$!!exty+LV$!!exty
                 set L$!!extz  = LR$!!extz+LV$!!extz
                 set E$!!extx  = ER$!!extx+EV$!!extx
                 set E$!!exty  = ER$!!exty+EV$!!exty
                 set E$!!extz  = ER$!!extz+EV$!!extz
                }
		if($filetype==$FIELDLINETYPE){\
		       #
		       }
                #
                if($filetype==$DTYPE){\
                 #rot180del r newr
                 #rot180delv vx newvx
                 #rot180delv vy newvy
		 # Accelerations
		 #accelerations 3 2
		}
                #
readb2pole 1    #
		define LOOPTYPE (1)
		rd 0_ener.dat
		da 0_ener.dat
		lines 4 10000000
		read {divbmax 114 divbavg 115 b2pole 116 polevolume 117}
		rd 0_loss.dat
		der t min1 td min1d
		set R1=b2pole/min1d
		define x2label "(B^2 r_h^2 c)/(\dot{M} c^2)"
		if($1==0){ ctype default pl 0 t R1 } else { ctype red plo 0 t R1 }
		#
readenerfull 10 #
		define filename 0_ener.dat
		if($?1==1){\
		 rd $filename $1
		 define ext $1
		 define extx $1x
		 define exty $1y
		 define extz $1z
		}\
		else{\
		 rd $filename
		 define exttype (0)
                 define ext " "
                 define extx x
                 define exty y
                 define extz z
                }
		da $filename
		lines 4 1000000
		read {nstep$!!ext 114 divbmax$!!ext 115 divbavg$!!ext 116 b2pole$!!ext 117 polevolume$!!ext 118}
		# for new collaboration
		#read {nstep$!!ext 114 divbmax$!!ext 115 divbavg$!!ext 116 b2pole$!!ext 117 polevolume$!!ext 118 iplus$!!ext 119 iminus$!!ext 120 volplus$!!ext 121 volminus$!!ext  122 massplus$!!ext 123 massminus$!!ext 124}
		#read {nstep$!!ext 114 divbmax$!!ext 115 divbavg$!!ext 116}
		#read {divbmax$!!ext 114 divbavg$!!ext 115}
		#read {nstep$!!ext 114 divbmax$!!ext 115 divbavg$!!ext 116}
		#read {nstep$!!ext 114 divbmax$!!ext 115 divbavg$!!ext 116 v2$!!extx 117 v2$!!exty 118 v2$!!extz 119 b2$!!extx 120 b2$!!exty 121 b2$!!extz 122}
		#read {nstep$!!ext 114 divbmax$!!ext 115 divbavg$!!ext 116 v2$!!extx 117 v2$!!exty 118 v2$!!extz 119 b2$!!extx 120 b2$!!exty 121 b2$!!extz 122 rp 123 enp 124 vpx 125 vpy 126 vpz 127 bpx 128 bpy 129 bpz 130}
