plc    17	# plc <file> <function> <type of plot=100,000,overlay=010,000,limits=000,001> <0 0 0 0>
		# this is a generic setup below
                if($?3 == 1) { define tobebits ($3) } else { define tobebits (0x0) }
                #defaults
                #expand 1.3
                #location 3500 17250 3500 31000
		#location 3500 31000 3500 31000
                #window 2 1 1 1
                myrd $1
                if('$tobebits'=='000') {
                 set thebits=0x0
		}
		if('$tobebits'!='000') {
                 set thebits=0x$tobebits
		}
		if(thebits & 0x001) {
                  shrink3 $2 x12 x22 $4 $5 $6 $7
		  set reallyx=x12new
		  set reallyy=x22new
		  set newfun=$2new
		 if($2=='r'){
		  set newfun=($2new)}
	
		   define xl (reallyx[0])
	   	   define xh (reallyx[$rnx-1])
		   define yl (reallyy[0])
		   define yh (reallyy[$rnx*$rny-1])
                   define txl ($xl)
                   define txh ($xh)
                   define tyl ($yl)
                   define tyh ($yh)
                }
		  if($PLANE==3) {
                   #image ($nx,$ny) $xl $xh $yl $yh
                   #image ($nx,$ny) $rxl $rxh $ryl $ryh
                   #image ($nx,$ny)
		   define realdx ($dx)
		   define realdy ($dy)
		   #
		   if(!(thebits & 0x001)) {
		   define xl ($Sx)
	   	   define xh ($Sx+$Lx)
		   define yl ($Sy)
		   define yh ($Sy+$Ly)
                   define txl ($xl)
                   define txh ($xh)
                   define tyl ($yl)
                   define tyh ($yh)
                   define rxl (0)
                   define rxh ($nx)
                   define ryl (0)
                   define ryh ($ny)
                   define rnx ($nx)
                   define rny ($ny)
		      #
                   #set reallyx=x1
                   #set reallyy=x2
                   #set reallyx=x12 # don't need for plc
                   #set reallyy=x22 # don't need for plc
		   #define WHICHLEV ($nz/2)
                   echo "pdimen here"
                   pdimen $2
                   set newfun=$2  if(k==$WHICHLEV) # for example
                 if($2==r){
                  set newfun=($2)}

		   set reallyx=x12 if(k==$WHICHLEV)
                   set reallyy=x22 if(k==$WHICHLEV)
		}
		   #
		 }
		 if($PLANE==2) {
		   #define WHICHLEV ($ny/2)
		   define realdx ($dx)
		   define realdy ($dz)
	   	   #
		   if(!(thebits & 0x001)) {
                   define xl ($Sx)
	   	   define xh ($Sx+$Lx)
		   define yl ($Sz)
		   define yh ($Sz+$Lz)
                   define txl ($xl)
                   define txh ($xh)
                   define tyl ($yl)
                   define tyh ($yh)
                   set newfun=$2 if(j==$WHICHLEV)
                 if($2=='r'){
                  set newfun=($2)}

                   define rxl (0)
                   define rxh ($nx)
                   define ryl (0)
                   define ryh ($nz)
                   define rnx ($nx)
                   define rny ($nz)
		   set reallyx=x12 if(j==$WHICHLEV)
                   set reallyy=x32 if(j==$WHICHLEV)
		}
		 }
		 if($PLANE==1) {
		   #define WHICHLEV ($nx/2)
		   define realdx ($dy)
		   define realdy ($dz)
		   if(!(thebits & 0x001)) {
                   define xl ($Sy)
	   	   define xh ($Sy+$Ly)
		   define yl ($Sz)
		   define yh ($Sz+$Lz)
                   define txl ($xl)
                   define txh ($xh)
                   define tyl ($yl)
                   define tyh ($yh)
                   define rxl (0)
                   define rxh ($ny)
                   define ryl (0)
                   define ryh ($nz)
                   define rnx ($ny)
                   define rny ($nz)
                   set newfun=$2 if(i==$WHICHLEV)   
                 if($2=='r'){
                  set newfun=($2)}

		   set reallyx=x22 if(i==$WHICHLEV)   
                   set reallyy=x32 if(i==$WHICHLEV)
		}
	   	   #
		 }
		#
		echo "got here"
		set ii=1,$rnx*$rny
		if($SKIPFACTOR>1) {
		 if( $SKIPFACTOR>$rnx/2-1 ||  $SKIPFACTOR>$rny/2-1) {
		  if($rny<$rnx) {
		   set temptemp=int($rny/2-1)
		   define SKIPFACTOR (temptemp)
		  }
		  if($rny>=$rnx) {
		   set temptemp=int($rnx/2-1)
		   define SKIPFACTOR (temptemp)
		  }
		 }
		 set ii=0,$rnx*$rny-1
		 set iix = ii%$rnx+1
		 set iiy = int(ii/$rnx)+1
		 set use = (int(iiy/$SKIPFACTOR) - iiy/$SKIPFACTOR == 0 && int(iix/$SKIPFACTOR) - iix/$SKIPFACTOR == 0  ) ? 1 : 0
		 set temptempnx=int($rnx/$SKIPFACTOR)
		 set temptempny=int($rny/$SKIPFACTOR)
		 define rnx (temptempnx)
		 define rny (temptempny)
		 set ii=0,$rnx*$rny-1
		 set ix = ii%$rnx
		 set iy = int(ii/$rnx)
		}
		if($SKIPFACTOR==1) {
		 set x=(ii - $rnx*int((ii-0.5)/$rnx) - 0.5)*$realdx+$txl
		 set y=((int((ii-0.5)/$rnx) + 0.5)*$realdy)+$tyl
		 set ix = int((x -$txl)/$realdx)
		 set iy = int((y -$tyl)/$realdy)
		 set use=ii/ii*1		 
		}
		#
		set newfuntemp=newfun if(use)
		set newfun=newfuntemp
		#set realx1=reallyx if(use)  # don't need for plc
		#set realy1=reallyy if(use) # don't need for plc
		#
		if(!(thebits & 0x010)) {
                 limits $txl $txh $tyl $tyh
		}
		image ($rnx,$rny) $txl $txh $tyl $tyh
		#
                set image[ix,iy] = newfun
                #limits $rxl $rxh $ryl $ryh
                #device ppm file1.ppm
                # default(0) to not overlay
                if(thebits & 0x010) {
                   define temptemptemp (1)
                }
		if(!(thebits & 0x010)) {
		 if($PLOTERASE) {
                #  erase
		 }
#		 ctype default
                 labeltime
                }
		# the above portion is a generic setup
		if($CONSTLEVELS==0) {
		 minmax min max
		}
		if($CONSTLEVELS!=0) {
		 if(constlevelshit==0) {
		  set constlevelshit=1
		  minmax min max
		 }
		}
		if($loginterprho==1) {
		  define min -3.7
		  define max -2
		}
                echo "min:"$min "max:"$max
                # define cres (15) # see plc
		if(($SOLIDCONTOURS==0)&&($min*$max < 0.)) {
                        echo "pt:1"
		        if($cres>0) {define delta (($max-$min)/$cres) }
		        if($cres<0) {define delta (-$cres) }
                        echo "delta:"$delta
                        if($min>-$delta) {
                          define min -$delta
                        }
			set lev=$min,-$delta,$delta
			levels lev
			ltype 0
#		        ctype $NEGCONTCOLOR
			contour		        
			#
                        if($max<$delta) {
                          define max $delta
                        }
			set lev=$delta,$max,$delta
			levels lev
			ltype 0
#		        ctype $POSCONTCOLOR
			contour
#		        ctype default
		}
		if(($min<=0)&&($max<=0)) {
                        echo "pt:2"
		        if($cres>0) { set lev=$min,$max,($max-$min)/$cres define delta (($max-$min)/$cres) }
		        if($cres<0) { set lev=$min,$max,-$cres define delta (-$cres) }
			levels lev
			#ltype 2
                        ltype 0
#		        ctype $NEGCONTCOLOR
			contour
		}
		if(($SOLIDCONTOURS==1)||($min>=0)&&($max>=0)) {
                        echo "pt:3"
		        if($cres>0) { set lev=$min,$max,($max-$min)/$cres define delta (($max-$min)/$cres) }
		        if($cres<0) { set lev=$min,$max,-$cres define delta (-$cres) }
			levels lev
			ltype 0
#		        ctype $POSCONTCOLOR
			contour
		}
		if($jet==1) {
		 set lev=.90,1.0,.01
		 levels lev
		 ctype yellow
		 ltype 0
		 contour
		 set lev=1.0,10.0,.01
		 levels lev
		 ctype blue
		 ltype 0
		 contour
		 ctype default
		}
                if(thebits & 0x010) {
                  define temptemptemp (0)
                }
		if(!(thebits & 0x010)) {
                 if(thebits & 0x100) {
                  ticksize -1 0 0 0
		  box
                  relocate (15000 31500)
                  label RADIAL SCALE NOT EXACTLY CORRECT
		 }
		 if(!(thebits & 0x100)) {
                  mybox2d
		  ticksize 0 0 0 0
		  limits $txl $txh $tyl $tyh		  
		 }
		 prepaxes x1 x2 $2
                 labelaxes 0
		}
                #
                #device X11
plcnoplot    17	# plc <file> <function> <type of plot=100,000,overlay=010,000,limits=000,001> <0 0 0 0>
		# this is a generic setup below
                if($?3 == 1) { define tobebits ($3) } else { define tobebits (0x0) }
                #defaults
                #expand 1.3
                #location 3500 17250 3500 31000
		#location 3500 31000 3500 31000
                #window 2 1 1 1
                myrd $1
                if('$tobebits'=='000') {
                 set thebits=0x0
		}
		if('$tobebits'!='000') {
                 set thebits=0x$tobebits
		}
		if(thebits & 0x001) {
                  shrink3 $2 x12 x22 $4 $5 $6 $7
		  set reallyx=x12new
		  set reallyy=x22new
		  set newfun=$2new
		   define xl (reallyx[0])
	   	   define xh (reallyx[$rnx-1])
		   define yl (reallyy[0])
		   define yh (reallyy[$rnx*$rny-1])
                   define txl ($xl)
                   define txh ($xh)
                   define tyl ($yl)
                   define tyh ($yh)
                }
		  if($PLANE==3) {
                   #image ($nx,$ny) $xl $xh $yl $yh
                   #image ($nx,$ny) $rxl $rxh $ryl $ryh
                   #image ($nx,$ny)
		   define realdx ($dx)
		   define realdy ($dy)
		   #
		   if(!(thebits & 0x001)) {
		   define xl ($Sx)
	   	   define xh ($Sx+$Lx)
		   define yl ($Sy)
		   define yh ($Sy+$Ly)
                   define txl ($xl)
                   define txh ($xh)
                   define tyl ($yl)
                   define tyh ($yh)
                   define rxl (0)
                   define rxh ($nx)
                   define ryl (0)
                   define ryh ($ny)
                   define rnx ($nx)
                   define rny ($ny)
		      #
                   #set reallyx=x1
                   #set reallyy=x2
                   #set reallyx=x12 # don't need for plc
                   #set reallyy=x22 # don't need for plc
		   #define WHICHLEV ($nz/2)
                   echo "pdimen here"
                   pdimen $2
                   set newfun=$2  if(k==$WHICHLEV) # for example
		   set reallyx=x12 if(k==$WHICHLEV)
                   set reallyy=x22 if(k==$WHICHLEV)
		}
		   #
		 }
		 if($PLANE==2) {
		   #define WHICHLEV ($ny/2)
		   define realdx ($dx)
		   define realdy ($dz)
	   	   #
		   if(!(thebits & 0x001)) {
                   define xl ($Sx)
	   	   define xh ($Sx+$Lx)
		   define yl ($Sz)
		   define yh ($Sz+$Lz)
                   define txl ($xl)
                   define txh ($xh)
                   define tyl ($yl)
                   define tyh ($yh)
                   set newfun=$2 if(j==$WHICHLEV)
                   define rxl (0)
                   define rxh ($nx)
                   define ryl (0)
                   define ryh ($nz)
                   define rnx ($nx)
                   define rny ($nz)
		   set reallyx=x12 if(j==$WHICHLEV)
                   set reallyy=x32 if(j==$WHICHLEV)
		}
		 }
		 if($PLANE==1) {
		   #define WHICHLEV ($nx/2)
		   define realdx ($dy)
		   define realdy ($dz)
		   if(!(thebits & 0x001)) {
                   define xl ($Sy)
	   	   define xh ($Sy+$Ly)
		   define yl ($Sz)
		   define yh ($Sz+$Lz)
                   define txl ($xl)
                   define txh ($xh)
                   define tyl ($yl)
                   define tyh ($yh)
                   define rxl (0)
                   define rxh ($ny)
                   define ryl (0)
                   define ryh ($nz)
                   define rnx ($ny)
                   define rny ($nz)
                   set newfun=$2 if(i==$WHICHLEV)   
		   set reallyx=x22 if(i==$WHICHLEV)   
                   set reallyy=x32 if(i==$WHICHLEV)
		}
	   	   #
		 }
		#
		echo "got here"
		set ii=1,$rnx*$rny
		if($SKIPFACTOR>1) {
		 if( $SKIPFACTOR>$rnx/2-1 ||  $SKIPFACTOR>$rny/2-1) {
		  if($rny<$rnx) {
		   set temptemp=int($rny/2-1)
		   define SKIPFACTOR (temptemp)
		  }
		  if($rny>=$rnx) {
		   set temptemp=int($rnx/2-1)
		   define SKIPFACTOR (temptemp)
		  }
		 }
		 set ii=0,$rnx*$rny-1
		 set iix = ii%$rnx+1
		 set iiy = int(ii/$rnx)+1
		 set use = (int(iiy/$SKIPFACTOR) - iiy/$SKIPFACTOR == 0 && int(iix/$SKIPFACTOR) - iix/$SKIPFACTOR == 0  ) ? 1 : 0
		 set temptempnx=int($rnx/$SKIPFACTOR)
		 set temptempny=int($rny/$SKIPFACTOR)
		 define rnx (temptempnx)
		 define rny (temptempny)
		 set ii=0,$rnx*$rny-1
		 set ix = ii%$rnx
		 set iy = int(ii/$rnx)
		}
		if($SKIPFACTOR==1) {
		 set x=(ii - $rnx*int((ii-0.5)/$rnx) - 0.5)*$realdx+$txl
		 set y=((int((ii-0.5)/$rnx) + 0.5)*$realdy)+$tyl
		 set ix = int((x -$txl)/$realdx)
		 set iy = int((y -$tyl)/$realdy)
		 set use=ii/ii*1		 
		}
		#
		set newfuntemp=newfun if(use)
		set newfun=newfuntemp
		#set realx1=reallyx if(use)  # don't need for plc
		#set realy1=reallyy if(use) # don't need for plc
		#
		if(!(thebits & 0x010)) {
                 limits $txl $txh $tyl $tyh
		}
		image ($rnx,$rny) $txl $txh $tyl $tyh
		#
                set image[ix,iy] = newfun
                #limits $rxl $rxh $ryl $ryh
                #device ppm file1.ppm
                # default(0) to not overlay
                if(thebits & 0x010) {
                   define temptemptemp (1)
                }
		if(!(thebits & 0x010)) {
		   if($PLOTERASE) {
                  #erase
		 }
		 #ctype default
                 #labeltime
                }
		# the above portion is a generic setup
		if($CONSTLEVELS==0) {
		 minmax min max
		}
		if($CONSTLEVELS!=0) {
		 if(constlevelshit==0) {
		  set constlevelshit=1
		  minmax min max
		 }
		}
		if($loginterprho==1) {
		  define min -3.7
		  define max -2
		}
                echo "min:"$min "max:"$max
                # define cres (15) # see plc
		if(($SOLIDCONTOURS==0)&&($min*$max < 0.)) {
                        echo "pt:1"
		        if($cres>0) {define delta (($max-$min)/$cres) }
		        if($cres<0) {define delta (-$cres) }
                        echo "delta:"$delta
                        if($min>-$delta) {
                          define min -$delta
                        }
			set lev=$min,-$delta,$delta
			levels lev
			ltype 0
		        ctype $NEGCONTCOLOR
			#contour		        
			#
                        if($max<$delta) {
                          define max $delta
                        }
			set lev=$delta,$max,$delta
			levels lev
			ltype 0
		        ctype $POSCONTCOLOR
			#contour
		        ctype default
		}
		if(($min<=0)&&($max<=0)) {
                        echo "pt:2"
		        if($cres>0) { set lev=$min,$max,($max-$min)/$cres define delta (($max-$min)/$cres) }
		        if($cres<0) { set lev=$min,$max,-$cres define delta (-$cres) }
			levels lev
			#ltype 2
                        ltype 0
		        ctype $NEGCONTCOLOR
			#contour
		}
		if(($SOLIDCONTOURS==1)||($min>=0)&&($max>=0)) {
                        echo "pt:3"
		        if($cres>0) { set lev=$min,$max,($max-$min)/$cres define delta (($max-$min)/$cres) }
		        if($cres<0) { set lev=$min,$max,-$cres define delta (-$cres) }
			levels lev
			ltype 0
		        ctype $POSCONTCOLOR
			contour
		}
		if($jet==1) {
		 set lev=.90,1.0,.01
		 levels lev
		 ctype yellow
		 ltype 0
		 contour
		 set lev=1.0,10.0,.01
		 levels lev
		 ctype blue
		 ltype 0
		 contour
		 ctype default
		}
                if(thebits & 0x010) {
                  define temptemptemp (0)
                }
		if(!(thebits & 0x010)) {
                 if(thebits & 0x100) {
                  ticksize -1 0 0 0
		  #box
                  #relocate (15000 31500)
                  #label RADIAL SCALE NOT EXACTLY CORRECT
		 }
		 if(!(thebits & 0x100)) {
                  #mybox2d
		  ticksize 0 0 0 0
		  limits $txl $txh $tyl $tyh		  
		 }
		 #prepaxes x1 x2 $2
                 #labelaxes 0
		}
                #
                #device X11
pls    17	# pls <file> <function> <type of plot=100,000,overlay=010,000,limits=000,001> <0 0 0 0>
                prepaxes x1 x2 $2
                if($?3 == 1) { define tobebits ($3) } else { define tobebits (0x0) }
                #defaults
                #expand 1.3
                #location 3500 17250 3500 31000
                location 3500 31000 3500 31000
                #window 2 1 1 1
                myrd $1
                if('$tobebits'=='000'){\
                 set thebits=0x0
		}\
                else{\
                 set thebits=0x$tobebits
		}
                if(thebits & 0x001){\
		  shrink3 $2 x12 x22 $4 $5 $6 $7
		  set newfun=$2new
		   if(thebits & 0x100){\
		         set inew=i if( (ixold>=$rxl)&(ixold<=$rxh)&(iyold>=$ryl)&(iyold<=$ryh))
		         set jnew=j if( (ixold>=$rxl)&(ixold<=$rxh)&(iyold>=$ryl)&(iyold<=$ryh))
		          set reallyx=inew if(k==$WHICHLEV)
		          set reallyy=jnew if(k==$WHICHLEV)
		   }\
		   else{\
		    set reallyx=x12new if(k==$WHICHLEV)
		    set reallyy=x22new if(k==$WHICHLEV)
                   }
                }\
                else{\
		 if($PLANE==3){\
		   define xl ($Sx)
	   	   define xh ($Sx+$Lx)
		   define yl ($Sy)
		   define yh ($Sy+$Ly)
                   define txl ($xl)
                   define txh ($xh)
                   define tyl ($yl)
                   define tyh ($yh)
                   define rxl (0)
                   define rxh ($nx)
                   define ryl (0)
                   define ryh ($ny)
                   define rnx ($nx)
                   define rny ($ny)
                   #image ($nx,$ny) $xl $xh $yl $yh
                   #image ($nx,$ny) $rxl $rxh $ryl $ryh
                   #image ($nx,$ny)
                   #set reallyx=x1
                   #set reallyy=x2
                   #set reallyx=x12 # don't need for plc
                   #set reallyy=x22 # don't need for plc
		   #define WHICHLEV ($nz/2)
                   set newfun=$2  if(k==$WHICHLEV) # for example
		   define realdx (($txh-$txl)/$nx)
		   define realdy (($tyh-$tyl)/$ny)
		   # define realdx ($dx)
		   # define realdy ($dy)
		   if(thebits & 0x100){\
		          set reallyx=i if(k==$WHICHLEV)
		          set reallyy=j if(k==$WHICHLEV)
		   }\
		   else{\
		    set reallyx=x12 if(k==$WHICHLEV)
		    set reallyy=x22 if(k==$WHICHLEV)
                   }
		#
		 }
		 if($PLANE==2){\
                   define xl ($Sx)
	   	   define xh ($Sx+$Lx)
		   define yl ($Sz)
		   define yh ($Sz+$Lz)
                   define txl ($xl)
                   define txh ($xh)
                   define tyl ($yl)
                   define tyh ($yh)
                   define rxl (0)
                   define rxh ($nx)
                   define ryl (0)
                   define ryh ($nz)
                   define rnx ($nx)
                   define rny ($nz)
		   #define WHICHLEV ($ny/2)
                   set newfun=$2 if(j==$WHICHLEV)
		   define realdx ($dx)
		   define realdy ($dz)
	   	   #
		   set reallyx=x12 if(j==$WHICHLEV)
                   set reallyy=x32 if(j==$WHICHLEV)
		 }
		 if($PLANE==1){\
                   define xl ($Sy)
	   	   define xh ($Sy+$Ly)
		   define yl ($Sz)
		   define yh ($Sz+$Lz)
                   define txl ($xl)
                   define txh ($xh)
                   define tyl ($yl)
                   define tyh ($yh)
                   define rxl (0)
                   define rxh ($ny)
                   define ryl (0)
                   define ryh ($nz)
                   define rnx ($ny)
                   define rny ($nz)
		   #define WHICHLEV ($nx/2)
                   set newfun=$2 if(i==$WHICHLEV)   
		   define realdx ($dy)
		   define realdy ($dz)
		   set reallyx=x22 if(i==$WHICHLEV)   
                   set reallyy=x32 if(i==$WHICHLEV)   
	   	   #
		 }		 
                }
		set ii=1,$rnx*$rny
		if($SKIPFACTOR>1){\
		 if( $SKIPFACTOR>$rnx/2-1 ||  $SKIPFACTOR>$rny/2-1){\
		  if($rny<$rnx){\
		   set temptemp=int($rny/2-1)
		   define SKIPFACTOR (temptemp)
		  }\
		  else{\
		   set temptemp=int($rnx/2-1)
		   define SKIPFACTOR (temptemp)
		  }
		 }
		 set ii=0,$rnx*$rny-1
		 set iix = ii%$rnx+1
		 set iiy = int(ii/$rnx)+1
		 set use = (int(iiy/$SKIPFACTOR) - iiy/$SKIPFACTOR == 0 && int(iix/$SKIPFACTOR) - iix/$SKIPFACTOR == 0  ) ? 1 : 0
		 set temptempnx=int($rnx/$SKIPFACTOR)
		 set temptempny=int($rny/$SKIPFACTOR)
		 define rnx (temptempnx)
		 define rny (temptempny)
		 set ii=0,$rnx*$rny-1
		 set ix = ii%$rnx
		 set iy = int(ii/$rnx)
		}\
		else{\
		 set x=(ii - $rnx*int((ii-0.5)/$rnx) - 0.5)*$realdx+$txl
		 set y=((int((ii-0.5)/$rnx) + 0.5)*$realdy)+$tyl
		 set ix = int((x -$txl)/$realdx)
		 set iy = int((y -$tyl)/$realdy)
		 set use=ii/ii*1		 
		}
		#
		set newfuntemp=newfun if(use)
		set newfun=newfuntemp
		set realx1=reallyx if(use)
		set realy1=reallyy if(use)
		#
		# sm image command loses first row of data, so add a buffer of 0s
		# occurs for positive large values in the 0->$rnx for tj=0, for whatever crappy buggy reason!
		#set truenewfun=1,$rnx*($rny+1),1
		#set truenewfun=0*truenewfun
		#set truerealx1=0*truenewfun
		#set truerealy1=0*truenewfun
		#do kk=$rny,$rnx*$rny-1,1{
		#   set truenewfun[$kk]=newfun[$kk-$rny]
		#   set truerealx1[$kk]=realx1[$kk-$rny]
		#   set truerealy1[$kk]=realy1[$kk-$rny]
		#}
		#define rny ($rny+$rnx)
		#set newfun=truenewfun
		#set realx1=truerealx1
		#set realx2=truerealx2
		#
		#
		#
		if(!(thebits & 0x010)){\
                 limits $txl $txh $tyl $tyh
		}
		#image ($nx,$ny) $xl $xh $yl $yh
		#image ($rnx,$rny)
		image ($rnx,$rny)
		#
                set image[ix,iy] = newfun
		#set image[ix,iy] = $2[ii-1]
                set realx=0,$rnx-1
                set realy=0,$rny-1
                do kk=0,$rnx-1,1{
                  set realx[$kk] = realx1[$kk]
                }
                do kk=0,$rny-1,1{
                  set realy[$kk] = realy1[$kk*$rnx]
                }
		#
                #
		if(!(thebits & 0x010)){\
                 limits 0 $rnx 0 $rny
		}
                #limits $rxl $rxh $ryl $ryh
                minmax min max echo $min $max
		if(!(thebits & 0x010)){\
                 limits 0 $rnx 0 $rny
		}
                #limits $rxl $rxh $ryl $ryh
                #device ppm file1.ppm
                # default(0) to not overlay
                if(thebits & 0x010){\
                   define temptemptemp (1)
                }\
                else{
                 erase
		 ctype default
                 labeltime
                }
                # else overlay
		load surfaces
		Viewpoint 50 50 0
		#Viewpoint 80 80 0
                #Viewpoint 20 20 0
		#Viewpoint 50 -130 0
                #surface 13 $min $max realx realy
                #Surface 13 $min $max realx realy
                #Surface 13 $min $max ix iy
                # default to not be log-log
		Surface 103 $min $max realx realy
		#
                define Lo_x ($txl)
                define Hi_x ($txh)
                define Lo_y ($tyl)
                define Hi_y ($tyh)
                box3
                labelaxes3 x1 x2 $2
                #device X11
		#
		#
		#
vpl	19	# eg. vpl dump0001.dat v  1 12 010  <0 0 0 0>
		# vpl <file/0> <vec>   <v_leng> <12=xy,13=xz> <log&&overlay&&limits> <limits>
		# log=1 0=normal
		# overlay=0, no overlay, 1=overlay
		# limits=0 normal limits, 1=specified at line <limits>
                # this is a generic setup below(modified for 2 extra head entries
		# can control how many vectors reduced by defining SKIPFACTOR
		# can control length of vector as uniform or not with UNITVECTOR=1/0
                if($?5 == 1) { define tobebits ($5) } else { define tobebits (0x0) }
                #defaults
                #expand 1.3
                #location 3500 17250 3500 31000
                #location 3500 31000 3500 31000
                #window 2 1 1 1
                myrd $1
                if('$tobebits'=='000'){\
                 set thebits=0x0
		}\
                else{\
                 set thebits=0x$tobebits
		}
                if(thebits & 0x001){\
		  shrink3 $2x x12 x22 $6 $7 $8 $9
		  if($4==12){\
		   shrink3 $2y x12 x22 $6 $7 $8 $9
		  }
		  if($4==13){\
		   shrink3 $2z x12 x22 $6 $7 $8 $9
		  }
		  set reallyx=x12new
		  set reallyy=x22new
		  set newfunx=$2xnew
		  set newfuny=$2ynew
		  set newfunz=$2znew
		  #
                  # do not change limits if overlay
		  if(!(thebits & 0x010)){\
                   limits $txl $txh $tyl $tyh
		  }
                  #image ($nx,$ny) $xl $xh $yl $yh
                  #image ($rnx,$rny)
                }\
                else{\
                   define xl ($Sx)
	   	   define xh ($Sx+$Lx)
		   define yl ($Sy)
		   define yh ($Sy+$Ly)
                   # do not change limits if overlay
		   if(!(thebits & 0x010)){\
                    limits $xl $xh $yl $yh
		   }
                   define txl ($xl)
                   define txh ($xh)
                   define tyl ($yl)
                   define tyh ($yh)
                   define rxl (0)
                   define rxh ($nx)
                   define ryl (0)
                   define ryh ($ny)
                   define rnx ($nx)
                   define rny ($ny)
                   #image ($nx,$ny) $xl $xh $yl $yh
                   #image ($nx,$ny) $rxl $rxh $ryl $ryh
                   #image ($nx,$ny)
		   define dx (($txh-$txl)/$nx)
		   define dy (($tyh-$tyl)/$ny)
		   set ii=1,$nx*$ny
                   set x=(ii - $nx*int((ii-0.5)/$nx) - 0.5)*$dx+$xl
		   set y=((int((ii-0.5)/$nx) + 0.5)*$dy)+$yl
		   set ix = int((x -$xl)/$dx)
		   set iy = int((y -$yl)/$dy)
                   #set reallyx=x1
                   #set reallyy=x2
                   set reallyx=x12
                   set reallyy=x22
                   set newfunx=$2x
		   if($4==12){\
		    set newfuny=$2y
	 	   }
		   if($4==13){\
		    set newfunz=$2z
	 	   }
	   	   # 
                }
                #limits $rxl $rxh $ryl $ryh
                #device ppm file1.ppm
                # default(0) to not overlay
                if(thebits & 0x010){\
                   define temptemptemp (1)
                }\
                else{
                 erase
		 ctype default
                 labeltime
                }
		#
		set ii=1,$rnx*$rny        
		if($SKIPFACTOR>1){\
		 #
		 set iy = int(ii/$rnx) + 1
		 set ix = ii - $rnx*(iy - 1)
		 set use = (int(iy/$SKIPFACTOR) - iy/$SKIPFACTOR == 0 && int(ix/$SKIPFACTOR) - ix/$SKIPFACTOR == 0) ? 1 : 0
		}\
		else{\
		        set use=ii/ii*1
		}
		set realx=reallyx if(use)
		set realy=reallyy if(use)
                if($4==13){\
			    set VVx = newfunx if(use)
		            set VVy = newfunz if(use)
			    }
		if($4==12){\
		       set VVx = newfunx if(use)
		       set VVy = newfuny if(use)
		       }
		set ang=180.*atan2(VVy,VVx)/pi
                # for cart coords only
                # spherical coord mapping not done here but in code using interpolation
		if(!($UNITVECTOR)){\
		 set leng=$3*sqrt(VVx**2 + VVy**2)
		}\
		else{\
		 set leng=$3
		}
                if(thebits & 0x100){\
                 ticksize 0 0 0 0
		 #limits (LG(realx)) realy
		 #limits (LG(realx)) realy
		   if(!(thebits & 0x010)){\
		        limits x1 x2
		     }
                 if(thebits & 0x010){\
		  define temptemp (0)
		 }\
		 else{\
		  box
		 }
		 #vfield (LG(realx)) realy leng ang
		 #vfield (LG(realx)) realy leng ang
		 #vfield x1 x2 leng ang
		 vfield x y leng ang
		 #set myx1=x1*1E3
		 #set myx2=x2*pi
		 #vfield myx1 myx2 leng ang
                }\
                else{\
                 ticksize 0 0 0 0
                 #vfield x1 x2 leng ang
                 vfield realx realy leng ang
                 if(thebits & 0x010){\
		  define temptemp (0)
		 }\
                 else{\
                  #mybox2d
                  box
 		  labelaxes 0
                  prepaxes x1 x2 $2
		 }
                }
                #
pl	18	# pl <file> <dir> <function> <logx=1000,0000,logy=0100,0000,overlay=0010,0000,limits=0000,0001> <0 0 0 0>
                prepaxes $2 $3
		lweight 4
                if($?4 == 1) { define tobebits ($4) } else { define tobebits (0x0) }
                #defaults
                myrd $1
                if( ('$tobebits'=='0000')||('$tobebits'=='0x0') ){\
                 set thebits=0x0
		}\
                else{\
                 set thebits=0x$tobebits
		}
		#
		#
                if(thebits & 0x1000){\
                 set rlx=LG(ABS($2))
		 if(thebits & 0x0001){\
		  set _newxl=LG($5)
		  define newxl (_newxl)
		  set _newxh=LG($6)
		  define newxh (_newxh)
		 }
                }\
                else{\
		 set rlx=$2
                 if(thebits & 0x0001){\
		  set _newxl=$5
		  define newxl (_newxl)
		  set _newxh=$6
		  define newxh (_newxh)
		 }
                }
                if(thebits & 0x0100){\
		 set rly=LG(ABS($3))
                 if(thebits & 0x0001){\
		  set _newyl=LG($7)
		  define newyl (_newyl)
		  set _newyh=LG($8)
		  define newyh (_newyh)
		 }
                }\
                else{\
		 set rly=$3
                 if(thebits & 0x0001){\
		  set _newyl=$7
		  define newyl (_newyl)
		  set _newyh=$8
		  define newyh (_newyh)
		 }
                }
	        if(!(thebits&0x0010)){\
		 if((thebits&0x1000)&&(thebits&0x0100)){\
		  #ticksize -.1 1 -1 10
		  ticksize -1 10 -1 10
		 }
		 if((thebits&0x1000)&&(!(thebits&0x0100))){\
		  ticksize -1 0 0 0
		 }
		 if((!(thebits&0x1000))&&(thebits&0x0100)){\
		  ticksize 0 0 -1 0
		 }
		 if((!(thebits&0x1000))&&(!(thebits&0x0100))){\
		  ticksize 0 0 0 0
		 }
		}
                # (ones->limits)
                # if ones==1 use input limits
                # else if ones==0 use standard limits
                # if tens=1 overlay
                if(thebits & 0x0010){\
		  # define your own ptype,ctype
                  define temptemptemp (0)
                }\
                else{\
                 if(thebits & 0x0001){\
		  limits  $newxl $newxh $newyl $newyh
                 }\
                 else{\
		  limits rlx rly
                 }
		 if($PLOTERASE){\
		  erase
		 }
		 labeltime
                }
		#
                mypoints rlx rly
                if(thebits & 0x0010){\
                  define temptemptemp (0)
                }\
                else{\
                  # use real names for labelaxes function
                  #labelaxes $2 $3
                  labelaxes 0
                  mybox1d 
                }
                # some non-general short cuts
plo     3       #
                pl $1 $2 $3 0010
pllim	7	#
                pl $1 $2 $3 0001 $4 $5 $6 $7
setlimits 6     # set limits for dotallslim
                  set lim2=0,5,1
                  set lim2[0]=$1
                  set lim2[1]=$2
                  set lim2[2]=$3
                  set lim2[3]=$4
                  set lim2[4]=$5
                  set lim2[5]=$6
                  define xl2 (lim2[0])
                  define xh2 (lim2[1])
                  define yl2 (lim2[2])
                  define yh2 (lim2[3])
                  define zl2 (lim2[4])
                  define zh2 (lim2[5])
                  #
gamgo           #
                cd /us2/jon/mhdt4/
                rdbasic 0 1 -1
                rd dump0026.dat          
                #rd dump0012.dat
                setlimits 0 6.28 3.14 3.14 0 1
		plflim 0 x1 p
                set it1=newfun
                setlimits 0 6.28 3.13 3.13 0 1
                plflim 0 x1 p
                set it2=newfun
                set itavg=0.5*(it1+it2)
                set myx=newx
                # gammie traces through y first!
                da /usr/u/gammie/rtwod/nonrel/dump010
                #da /usr/u/gammie/rtwod/nonrel/dump005
                lines 2 10000000
                read {x 1 y 2 rho 3 ux 4 uy 5 uz 6 Bx 7 By 8 Bz 9 u 10}
                define nx (512)
                define ny (512)
                set x12=pi+y # really y
                set x22=pi+x # really x
                set p=2.0/3.0*u
                setlimits 3.14 3.14 0 6.28 0 1
                ctype red plflim 0 x2 p
                set git1=newfun
                set gy=newx
                setlimits 3.13 3.13 0 6.28 0 1
                ctype red plflim 0 x2 p
                set git2=newfun
                set gavg=0.5*(git1+git2)
                #
                ptype 2 0
                ctype default pl 0 myx itavg
                points myx itavg
                ptype 4 0
                ctype red plo 0 gy gavg
                points gy gavg
plflim 17       #
                #plflim	17 : plflim <file> <dir> <x> <y> <f(x,y)> <flimit=0,1> <logx=100,000,logy=010,000,overlay=001,000>
                #eg. setlimits  1.0 1.4 0 3.14159 0 1 plflim dump0000.dat x1 x12 x22 r 0 100
                #eg. setlimits  0 4 1.55 1.58 0 1 plflim dump0000.dat x1 x12 x22 r 0 000
                # normal strong limits
                # setlimits  1 1.2 1.55 1.58 0 1 plflim dump0000.dat x1 x12 x22 r 0 000
                #
                if($?6 == 1) { define tobebits2 ($6) } else { define tobebits2 (0x0) }
                if($?7 == 1) { define tobebits ($7) } else { define tobebits (0x0) }
                if( ('$tobebits2'=='000')||('$tobebits2'=='0x0') ){\
                 set thebits2=0x0
		}\
                else{\
                 set thebits2=0x$tobebits2
		}
                if( ('$tobebits'=='000')||('$tobebits'=='0x0') ){\
                 set thebits=0x0
		}\
                else{\
                 set thebits=0x$tobebits
		}
                myrd $1
		shrink3 $5 $3 $4 $xl2 $xh2 $yl2 $yh2
		set reallyx=$3new
		set reallyy=$4new
		set newfun=$5new
		#
		if(('$2'=='x1')||('$2'=='x12')){\
                 set newx=reallyx
		 define newxl2 ($xl2)
		 define newxh2 ($xh2)
		}
		if(('$2'=='x2')||('$2'=='x22')){\
		 set newx=reallyy
		 define newxl2 ($yl2)
		 define newxh2 ($yh2)
		}
      	        if(thebits2 & 0x1){\
		 if(thebits & 0x111){
		  define newbits ($71)
                  pl 0 newx newfun $71 $newxl2 $newxh2 $zl2 $zh2
		 }\
                 else{\
                  pl 0 newx newfun 0001 $newxl2 $newxh2 $zl2 $zh2
		 }
		}\
                else{\
		 if(thebits & 0x111){
                  pl 0 newx newfun $70
		 }\
		 else{\
                  pl 0 newx newfun
		 }
		}
                #
pld	4	#
                #defaults
                myrd $1
                myrd $2 2
                set rdiff = r-r2
                set endiff = en-en2
                set sediff = se-se2
                set potdiff = pot-pot2
                set pdiff = p-p2
                set csdiff = cs-cs2
                set vxdiff = vx-v2x
                set vydiff = vy-v2y
                set vzdiff = vz-v2z
                set bxdiff = bx-b2x
                set bydiff = by-b2y
                set bzdiff = bz-b2z
                pl 0 $3 $4
pldo	4	#
                myrd $1
                myrd $2 2
                set rdiff = r-r2
                set endiff = en-en2
                set sediff = se-se2
                set potdiff = pot-pot2
                set pdiff = p-p2
                set csdiff = cs-cs2
                set vxdiff = vx-v2x
                set vydiff = vy-v2y
                set vzdiff = vz-v2z
		plo 0 $3 $4
pldlim	8	#
                #defaults
                myrd $1
                rd $2 2
                set rdiff = r-r2
                set endiff = en-en2
                set pdiff = p-p2
                set vxdiff = vx-v2x
                set vydiff = vy-v2y
                set vzdiff = vz-v2z
                set potdiff = pot-pot2
		limits $5 $6 $7 $8
		erase
                labeltime
                ctype default
		box
		ptype 4 0
		points $3 $4
		connect $3 $4
		xla $3
		yla $4
pldlimo	8	#
                myrd $1
                rd $2 2
                set rdiff = r-r2
                set endiff = en-en2
                set pdiff = p-p2
                set vxdiff = vx-vx2
                set vydiff = vy-vy2
                set vzdiff = vz-vz2
                set potdiff = pot-pot2
                ctype blue
		ptype 4 0
		points $3 $4
		connect $3 $4
		xla $3
		yla $4
pla	2	#
                #defaults
                myrd $1
		set ii=1,$nx*$ny
                set x=(ii - $nx*int((ii-0.5)/$nx) - 0.5)*$dx
		set y=(int((ii-0.5)/$nx) + 0.5)*$dy
		image ($nx,$ny) -.5 .5 -.5 .5
		set ix = int(x/$dx)
		set iy = int(y/$dy)
		set image (ix,iy) = $2[ii-1]
		#
		minmax min max echo $min $max
		#
		define delta (($max-$min)/10.)
		set lev=$min,-$delta,$delta
		levels lev
		ltype 2
		contour
		#
		set lev=$delta,$max,$delta
		levels lev
		ltype 0
		contour
		#
                # device postencap /us1/jon/compare/fiducial/r2d.eps
avgtimeinit  1  #
		set $1=1,$nx*$ny*$nz,1
		set $1=0*$1
avgtimesub   2  # avgtimesub rtime r
		set $1=$1+$2
avgtimefinish 1     #
		set $1=$1/numtotal
                #
avgtime   3	# avgtime (e.g. avgtime 'dump' start end)
		define DOCALCS 1
		define DOSTRESS 1
		define DOACC 0
		define FFV 0
		#defaults
		rdnumd
                set h1=$1
		            set h0np='np'
                set h3='.dat'
		#
		avgtimeinit rtime
		avgtimeinit entime
		avgtimeinit betime
		avgtimeinit cs2time
		avgtimeinit b2time
		avgtimeinit stime
		avgtimeinit vtimex
		avgtimeinit vtimey
		avgtimeinit vtimez
		if($fullvec){\
                       avgtimeinit btimex
                       avgtimeinit btimey
                       avgtimeinit btimez
                       avgtimeinit vatime
                       avgtimeinit vatimex
                       avgtimeinit FME1timex
                       avgtimeinit FME1timey
                       avgtimeinit FME1timez
                       avgtimeinit FME2timex
                       avgtimeinit FME2timey
                       avgtimeinit FME2timez
		    }   
		    avgtimeinit FPEtimex
		    avgtimeinit FPEtimey
		    avgtimeinit FPEtimez
		    avgtimeinit FEtimex
		    avgtimeinit FEtimey
		    avgtimeinit FEtimez
		if($DOSTRESS){\
                       avgtimeinit Fltime
                       avgtimeinit Flrtime
                       avgtimeinit Flmtime
                       avgtimeinit Flatime
                       avgtimeinit Flartime
                       avgtimeinit Flamtime
                       avgtimeinit Mftimex
                       avgtimeinit Mftimey
                       avgtimeinit Lftimex
                       avgtimeinit Lftimey
		    }
		if($npdone){\
                       avgtimeinit sig11time
                       avgtimeinit sig12time
                       avgtimeinit sig13time
                       avgtimeinit sig22time
                       avgtimeinit sig23time
                       avgtimeinit sig33time
                       avgtimeinit Fvzdxtime
                       avgtimeinit Fvzdytime
                       avgtimeinit Fmdxtime
                       avgtimeinit Fmdytime
		    }
		if($DOACC){\
		           avgtimeinit atvtimex
		           avgtimeinit atvtimey
		           avgtimeinit atvtimez
		           avgtimeinit aphitimex
		           avgtimeinit aphitimey
		           avgtimeinit aptimex
		           avgtimeinit aptimey
		           avgtimeinit abptimex
		           avgtimeinit abptimey
		           avgtimeinit atbtimex
		           avgtimeinit atbtimey
		           avgtimeinit atbtimez
		   }        
		 if($FFV){\
		           avgtimeinit avtimex
		           avgtimeinit avtimey
		           avgtimeinit avtimez
		           avgtimeinit aptimex
		           avgtimeinit aptimey
		           avgtimeinit aptimez
		           avgtimeinit agtimex
		           avgtimeinit agtimey
		           avgtimeinit agtimez
		           avgtimeinit acbtimex
		           avgtimeinit acbtimey
		           avgtimeinit acbtimez
		           avgtimeinit ab2timex
		           avgtimeinit ab2timey
		           avgtimeinit ab2timez
		           avgtimeinit abtimex
		           avgtimeinit abtimey
		           avgtimeinit abtimez
		           avgtimeinit fketimex
		           avgtimeinit fketimey
		           avgtimeinit fketimez
		           avgtimeinit fhtimex
		           avgtimeinit fhtimey
		           avgtimeinit fhtimez
		           avgtimeinit fgtimex
		           avgtimeinit fgtimey
		           avgtimeinit fgtimez
		           avgtimeinit fbtimex
		           avgtimeinit fbtimey
		           avgtimeinit fbtimez
		           avgtimeinit fb2timex
		           avgtimeinit fb2timey
		           avgtimeinit fb2timez
		           avgtimeinit fmtimex
		           avgtimeinit fmtimey
		           avgtimeinit fmtimez
		           avgtimeinit fptimexx
		           avgtimeinit fptimexy
		           avgtimeinit fptimexz
		           avgtimeinit fptimeyx
		           avgtimeinit fptimeyy
		           avgtimeinit fptimeyz
		           avgtimeinit fptimezx
		           avgtimeinit fptimezy
		           avgtimeinit fptimezz
		           avgtimeinit fbtimexx
		           avgtimeinit fbtimexy
		           avgtimeinit fbtimexz
		           avgtimeinit fbtimeyx
		           avgtimeinit fbtimeyy
		           avgtimeinit fbtimeyz
		           avgtimeinit fbtimezx
		           avgtimeinit fbtimezy
		           avgtimeinit fbtimezz
		           avgtimeinit vmtime
		           avgtimeinit vketimex
		           avgtimeinit vketimey
		           avgtimeinit vketimez
		           avgtimeinit vietime
		           avgtimeinit vgetime
		           avgtimeinit vbetimex
		           avgtimeinit vbetimey
		           avgtimeinit vbetimez
		           avgtimeinit vptimex
		           avgtimeinit vptimey
		           avgtimeinit vptimez
		           avgtimeinit vltimex
		           avgtimeinit vltimey
		           avgtimeinit vltimez
		           avgtimeinit vbtimex
		           avgtimeinit vbtimey
		           avgtimeinit vbtimez
		 }
                #
		#set numstart=20
		#set numend=50
		#set numstart=$NUMDUMPS/2
		#set numend=$NUMDUMPS-1
		#set numstart=2*$NUMDUMPS/3
		#set numend=3*$NUMDUMPS/5
		set numstart=$2
		set numend=$3
                set numtotal=numend-numstart+1
                do ii=numstart,numend,1 {
                  set h2=sprintf('%04d',$ii)
                  set _fname=h1+h2+h3
		              set _fnpname=h0np+h1+h2+h3
                  define filename (_fname)
                  define npfilename (_fnpname)
		  rd $filename
		  #
		  avgtimesub rtime r
		  avgtimesub entime en
		  avgtimesub betime Be
		  avgtimesub cs2time (cs*cs)
		  avgtimesub stime entropy
		  avgtimesub vtimex vx
		  avgtimesub vtimey vy
		  avgtimesub vtimez vz
		  avgtimesub FPEtimex FPEx
		  avgtimesub FPEtimey FPEy
		  avgtimesub FPEtimez FPEz
		  avgtimesub FEtimex FEx
		  avgtimesub FEtimey FEy
		  avgtimesub FEtimez FEz
		  if($fullvec){\
			 avgtimesub btimex bx
			 avgtimesub btimey by
			 avgtimesub btimex bz
			 avgtimesub vatime va
			 avgtimesub vatimex vax
			 avgtimesub b2time (bx*bx+by*by+bz*bz)
                       avgtimesub FME1timex FME1x
                       avgtimesub FME1timey FME1y
                       avgtimesub FME1timez FME1z
                       avgtimesub FME2timex FME2x
                       avgtimesub FME2timey FME2y
                       avgtimesub FME2timez FME2z
		      }
		  if($DOSTRESS){\
                         avgtimesub Mftimex Fmx
                         avgtimesub Mftimey Fmy
                         avgtimesub Lftimex Fvzx
                         avgtimesub Lftimey Fvzy
		                     avgtimesub Flrtime (r*vx*vz)
		                     avgtimesub Fltime (r*vx*vz-bz*bx)
			                   avgtimesub Flmtime (-bz*bx)		                     
		                     avgtimesub Flartime (x12*x12*x12*g42*r*vx*vz)
		                     avgtimesub Flatime (x12*x12*x12*g42*(r*vx*vz-bz*bx))
			                   avgtimesub Flamtime (x12*x12*x12*g42*(-bz*bx))
		      }
		  if($npdone){\
                         avgtimesub sig11time sig11
                         avgtimesub sig12time sig12
                         avgtimesub sig13time sig13
                         avgtimesub sig22time sig22
                         avgtimesub sig23time sig23
                         avgtimesub sig33time sig33
                         avgtimesub Fvzdtimex Fvzdx
                         avgtimesub Fvzdtimey Fvzdy
                         avgtimesub Fmdtimex Fmdx
                         avgtimesub Fmdtimey Fmdy
		      }
		      if($FFV){\
		           rd $npfilename
		           avgtimeinit avtimex
		           avgtimeinit avtimey
		           avgtimeinit avtimez
		           avgtimeinit aptimex
		           avgtimeinit aptimey
		           avgtimeinit aptimez
		           avgtimeinit agtimex
		           avgtimeinit agtimey
		           avgtimeinit agtimez
		           avgtimeinit acbtimex
		           avgtimeinit acbtimey
		           avgtimeinit acbtimez
		           avgtimeinit ab2timex
		           avgtimeinit ab2timey
		           avgtimeinit ab2timez
		           avgtimeinit abtimex
		           avgtimeinit abtimey
		           avgtimeinit abtimez
		           avgtimeinit fketimex
		           avgtimeinit fketimey
		           avgtimeinit fketimez
		           avgtimeinit fhtimex
		           avgtimeinit fhtimey
		           avgtimeinit fhtimez
		           avgtimeinit fgtimex
		           avgtimeinit fgtimey
		           avgtimeinit fgtimez
		           avgtimeinit fbtimex
		           avgtimeinit fbtimey
		           avgtimeinit fbtimez
		           avgtimeinit fb2timex
		           avgtimeinit fb2timey
		           avgtimeinit fb2timez
		           avgtimeinit fmtimex
		           avgtimeinit fmtimey
		           avgtimeinit fmtimez
		           avgtimeinit fptimexx
		           avgtimeinit fptimexy
		           avgtimeinit fptimexz
		           avgtimeinit fptimeyx
		           avgtimeinit fptimeyy
		           avgtimeinit fptimeyz
		           avgtimeinit fptimezx
		           avgtimeinit fptimezy
		           avgtimeinit fptimezz
		           avgtimeinit fbtimexx
		           avgtimeinit fbtimexy
		           avgtimeinit fbtimexz
		           avgtimeinit fbtimeyx
		           avgtimeinit fbtimeyy
		           avgtimeinit fbtimeyz
		           avgtimeinit fbtimezx
		           avgtimeinit fbtimezy
		           avgtimeinit fbtimezz
		           avgtimeinit vmtime
		           avgtimeinit vketimex
		           avgtimeinit vketimey
		           avgtimeinit vketimez
		           avgtimeinit vietime
		           avgtimeinit vgetime
		           avgtimeinit vbetimex
		           avgtimeinit vbetimey
		           avgtimeinit vbetimez
		           avgtimeinit vptimex
		           avgtimeinit vptimey
		           avgtimeinit vptimez
		           avgtimeinit vltimex
		           avgtimeinit vltimey
		           avgtimeinit vltimez
		           avgtimeinit vbtimex
		           avgtimeinit vbtimey
		           avgtimeinit vbtimez
		}
		if($DOACC){\
		           accelerations 3 2
		           avgtimesub atvtimex atvx
		           avgtimesub atvtimey atvy
		           avgtimesub atvtimez atvz
		           avgtimesub aphitimex aphix
		           avgtimesub aphitimey aphiy
		           avgtimesub aptimex apx
		           avgtimesub aptimey apy
		           avgtimesub abptimex abpx
		           avgtimesub abptimey abpy
		           avgtimesub atbtimex atbx
		           avgtimesub atbtimey atby
		           avgtimesub atbtimez atbz

		   }        
		}
	        #
		avgtimefinish rtime
		avgtimefinish entime
		avgtimefinish betime
		avgtimefinish cs2time
		avgtimefinish b2time
		avgtimefinish stime
		avgtimefinish vtimex
		avgtimefinish vtimey
		avgtimefinish vtimez
		avgtimefinish FPEtimex
		avgtimefinish FPEtimey
		avgtimefinish FPEtimez
		avgtimefinish FEtimex
		avgtimefinish FEtimey
		avgtimefinish FEtimez
		if($fullvec){\
                       avgtimefinish btimex
                       avgtimefinish btimey
                       avgtimefinish btimez
                       avgtimefinish vatime
                       avgtimefinish vatimex
                       avgtimefinish FME1timex
                       avgtimefinish FME1timey
                       avgtimefinish FME1timez
                       avgtimefinish FME2timex
                       avgtimefinish FME2timey
                       avgtimefinish FME2timez
		    }   
		if($DOSTRESS){\
                       avgtimefinish Fltime
                       avgtimefinish Flrtime
                       avgtimefinish Flmtime
                       avgtimefinish Flatime
                       avgtimefinish Flartime
                       avgtimefinish Flamtime
                       avgtimefinish Mftimex
                       avgtimefinish Mftimey
                       avgtimefinish Lftimex
                       avgtimefinish Lftimey
		    }
		if($npdone){\
                       avgtimefinish sig11time
                       avgtimefinish sig12time
                       avgtimefinish sig13time
                       avgtimefinish sig22time
                       avgtimefinish sig23time
                       avgtimefinish sig33time
                       avgtimefinish Fvzdxtime
                       avgtimefinish Fvzdytime
                       avgtimefinish Fmdxtime
                       avgtimefinish Fmdytime
		    }
		 if($DOACC){\
		           avgtimefinish atvtimex
		           avgtimefinish atvtimey
		           avgtimefinish atvtimez
		           avgtimefinish aphitimex
		           avgtimefinish aphitimey
		           avgtimefinish aptimex
		           avgtimefinish aptimey
		           avgtimefinish abptimex
		           avgtimefinish abptimey
		           avgtimefinish atbtimex
		           avgtimefinish atbtimey
		           avgtimefinish atbtimez
		 }
		      if($FFV){\
		           avgtimefinish avtimex
		           avgtimefinish avtimey
		           avgtimefinish avtimez
		           avgtimefinish aptimex
		           avgtimefinish aptimey
		           avgtimefinish aptimez
		           avgtimefinish agtimex
		           avgtimefinish agtimey
		           avgtimefinish agtimez
		           avgtimefinish acbtimex
		           avgtimefinish acbtimey
		           avgtimefinish acbtimez
		           avgtimefinish ab2timex
		           avgtimefinish ab2timey
		           avgtimefinish ab2timez
		           avgtimefinish abtimex
		           avgtimefinish abtimey
		           avgtimefinish abtimez
		           avgtimefinish fketimex
		           avgtimefinish fketimey
		           avgtimefinish fketimez
		           avgtimefinish fhtimex
		           avgtimefinish fhtimey
		           avgtimefinish fhtimez
		           avgtimefinish fgtimex
		           avgtimefinish fgtimey
		           avgtimefinish fgtimez
		           avgtimefinish fbtimex
		           avgtimefinish fbtimey
		           avgtimefinish fbtimez
		           avgtimefinish fb2timex
		           avgtimefinish fb2timey
		           avgtimefinish fb2timez
		           avgtimefinish fmtimex
		           avgtimefinish fmtimey
		           avgtimefinish fmtimez
		           avgtimefinish fptimexx
		           avgtimefinish fptimexy
		           avgtimefinish fptimexz
		           avgtimefinish fptimeyx
		           avgtimefinish fptimeyy
		           avgtimefinish fptimeyz
		           avgtimefinish fptimezx
		           avgtimefinish fptimezy
		           avgtimefinish fptimezz
		           avgtimefinish fbtimexx
		           avgtimefinish fbtimexy
		           avgtimefinish fbtimexz
		           avgtimefinish fbtimeyx
		           avgtimefinish fbtimeyy
		           avgtimefinish fbtimeyz
		           avgtimefinish fbtimezx
		           avgtimefinish fbtimezy
		           avgtimefinish fbtimezz
		           avgtimefinish vmtime
		           avgtimefinish vketimex
		           avgtimefinish vketimey
		           avgtimefinish vketimez
		           avgtimefinish vietime
		           avgtimefinish vgetime
		           avgtimefinish vbetimex
		           avgtimefinish vbetimey
		           avgtimefinish vbetimez
		           avgtimefinish vptimex
		           avgtimefinish vptimey
		           avgtimefinish vptimez
		           avgtimefinish vltimex
		           avgtimefinish vltimey
		           avgtimefinish vltimez
		           avgtimefinish vbtimex
		           avgtimefinish vbtimey
		           avgtimefinish vbtimez
		}
		# turn on calculations in case turned off
		define DOCALCS 1      
                #
animvpl	19	# animvpl 'dump' v 1 12 10  <0 0 0 0>
                #defaults
                if($?5 == 0) { define numsend (4) }\
                else{\
                  if($?6 == 1) { define numsend (6) } else { define numsend 5 }
                }
		rdnumd
                set h1=$1
                set h3='.dat'
                do ii=0,$NUMDUMPS-1,$ANIMSKIP {
                  set h2=sprintf('%04d',$ii)
                  set _fname=h1+h2+h3
                  define filename (_fname)
                  if($numsend==6){\
 		   vpl $filename $2 $3 $4 $5 $6 $7 $8 $9
                  }\
                  else{\
                   if($numsend==5){\
 		    vpl $filename $2 $3 $4 $5
                   }\
                   else{\
                    if($numsend==4){\
 		     vpl $filename $2 $3 $4
                    }
                   }
                  }
                  #delay loop
                  #set jj=0
                  #while {jj<1} {set jj=jj+1}
		}
animvplp 19	# animvplp 'dump' v 1 4 10  <0 0 0 0>
                #defaults
		rdnumd
                set h1=$1
                set h3='.dat'
                do ii=0,$NUMDUMPS-1,$ANIMSKIP {
                  set h2=sprintf('%04d',$ii)
                  set _fname=h1+h2+h3
                  define filename (_fname)
		  vplp $filename $2 $3 $4 $5 $6 $7 $8 $9
                  #delay loop
                  #set jj=0
                  #while {jj<1} {set jj=jj+1}
		}
		#
animplc	17	# animplc 'dump' r 000 <0 0 0 0>
                if($?3 == 0) { define numsend (2) }\
                else{\
                  if($?4 == 1) { define numsend (4) } else { define numsend 3 }
                }
                #defaults
		define PLANE (3)
		define WHICHLEV (0)
		#rdnumd
                set h1=$1
                set h3='.dat'
		set endanim=$NUMDUMPS-1
		#set endanim=$NUMFLDUMPS-1 
                do ii=0,endanim,$ANIMSKIP {
		  if($GAMMIE==0){ set h2=sprintf('%04d',$ii) set _fname=h1+h2+h3 }
		  if($GAMMIE==1){ set h2=sprintf('%04d',$ii) set _fname=h1+h2 }
                  
                  define filename (_fname)
                  if($numsend==2){ plc  $filename $2}\
                  else{\
                   if($numsend==3){  plc  $filename $2 $3}\
                   else{\
                    if($numsend==4){ plc  $filename $2 $3 $4 $5 $6 $7}
                   }
                  }
                  #delay loop
                  #set jj=0
                  #while {jj<1} {set jj=jj+1}
		}
		#
animplc2 18	# animplc 'dump' r r2 000 <0 0 0 0>
10.1.1.238
                if($?4 == 0) { define numsend (2) }\
                else{\
                  if($?5 == 1) { define numsend (4) } else { define numsend 3 }
                }
                #defaults
		define PLANE (3)
		define WHICHLEV (0)
		rdnumd
                set h1=$1
                set h3='.dat'
		set endanim=$NUMDUMPS-1
		#set endanim=$NUMFLDUMPS-1 
                do ii=0,endanim,$ANIMSKIP {
                  set h2=sprintf('%04d',$ii)
                  set _fname=h1+h2+h3
                  define filename (_fname)
		  plc  $filename $2 000
		  ctype blue
		  plc  0 $3 010
		  ctype default
                  #delay loop
                  #set jj=0
                  #while {jj<1} {set jj=jj+1}
		}
		#
animplcd 17	# animplc 'dump' r 000 <0 0 0 0>
                if($?3 == 0) { define numsend (2) }\
                else{\
                  if($?4 == 1) { define numsend (4) } else { define numsend 3 }
                }
                #defaults
		define PLANE (3)
		define WHICHLEV (0)
		rdnumd
                set h1=$1
                set h3='.dat'
		set endanim=$NUMDUMPS-1
		#set endanim=$NUMFLDUMPS-1
		#set floorlistmin=0,endanim,1
		#set floorlistmax=0,endanim,1
                do ii=0,endanim,$ANIMSKIP {
                  set h2=sprintf('%04d',$ii)
                  set _fname=h1+h2+h3
                  define filename (_fname)
                  if($numsend==2){ plc  $filename $2}\
                  else{\
                   if($numsend==3){  plc  $filename $2 $3}\
                   else{\
                    if($numsend==4){ plc  $filename $2 $3 $4 $5 $6 $7}
                   }
                  }
		   #set floorlistmin[$ii]=$min
		   #set floorlistmax[$ii]=$max
                  #delay loop
                  #set jj=0
                  #while {jj<1} {set jj=jj+1}
		}
		#
animzplc 17	# animzplc dump0000.dat r 000 <0 0 0 0>
                if($?3 == 0) { define numsend (2) }\
                else{\
                  if($?4 == 1) { define numsend (4) } else { define numsend 3 }
                }
                #defaults
		set startanim=0
		set endanim=$nz-3
		define PLANE 3
		myrd $1
		define ii (startanim)
		#do ii=startanim,endanim,$ANIMSKIP {\
		while{1==1}{\
		  set iitemp=$ii
		  print '%d\n' {iitemp}
		  define WHICHLEV ($ii)
                  if($numsend==2){ plc  0 $2}\
                  else{\
                   if($numsend==3){  plc  0 $2 $3}\
                   else{\
                    if($numsend==4){ plc  0 $2 $3 $4 $5 $6 $7}
                   }
                  }
		  set iitemp=$ii+1
		  define ii (iitemp)
		  if($ii==endanim){ define ii (startanim) }
                  #delay loop
                  set jj=0
                  while {jj<10000} {set jj=jj+1}
		}
		#
animzpls 17	# animzplc dump0000.dat r 000 <0 0 0 0>
                if($?3 == 0) { define numsend (2) }\
                else{\
                  if($?4 == 1) { define numsend (4) } else { define numsend 3 }
                }
                #defaults
		set startanim=0
		set endanim=$nz-3
		define PLANE 3
		myrd $1
		define ii (startanim)
		#do ii=startanim,endanim,$ANIMSKIP {\
		while{1==1}{\
		  set iitemp=$ii
		  define WHICHLEV ($ii)
		  print '%d\n' {iitemp}		  
                  if($numsend==2){ pls  0 $2}\
                  else{\
                   if($numsend==3){  pls  0 $2 $3}\
                   else{\
                    if($numsend==4){ pls  0 $2 $3 $4 $5 $6 $7}
                   }
                  }
		  set iitemp=$ii+1
		  define ii (iitemp)
		  if($ii==endanim){ define ii (startanim) }
                  #delay loop
                  set jj=0
                  while {jj<10000} {set jj=jj+1}
		}
		#
animplflc 17	# animplflc 'fldump' phibz 000 <0 0 0 0>
                if($?3 == 0) { define numsend (2) }\
                else{\
                  if($?4 == 1) { define numsend (4) } else { define numsend 3 }
                }
                #defaults
		rdnumd
                set h1=$1
                set h3='.dat'
		set h32='.ppm'
		define ANIMSKIP 1
		define PLOTERASE 0
		define CONSTLEVELS 0
		define SKIPFACTOR 1
		define cres 256
		set constlevelshit=0
		#set constlevelshit=1
		#define max 0
		#define min -3
		#define min -2.98509     # make center touch with cres=100
		set start=0   # don't restart unless min/max are actually the same(or set the same)
		#set finish=$NUMFLDUMPS-1
		set finish=$NUMDUMPS-1
                do ii=start,finish,$ANIMSKIP {
                  set h2=sprintf('%04d',$ii)
                  set _fname=h1+h2+h3
                  define filename (_fname)
		  set _fname2=h1+h2+h32
		  define filename2 (_fname2)
		  #device ppm ./fppm/$filename2
		  device ppm ./fppm/$filename2
		  if(1){\
                  if($numsend==2){ plc  $filename $2}\
                  else{\
                   if($numsend==3){  plc  $filename $2 $3}\
                   else{\
                    if($numsend==4){ plc  $filename $2 $3 $4 $5 $6 $7}
                   }
                  }
		  } else{\
		         echo pl
                  if($numsend==2){ pl  $filename x1 $2}\
                  else{\
                   if($numsend==3){  pl  $filename x1 $2 $3}\
                   else{\
                    if($numsend==4){ pl  $filename x1 $2 $3 $4 $5 $6 $7}
                   }
                  }
                  }
		  device X11
                  #delay loop
                  #set jj=0
                  #while {jj<1} {set jj=jj+1}
		}
		define PLOTERASE 1
		set constlevelshit=0
		define CONSTLEVELS 0
		#
animpls	17	#
                if($?3 == 0) { define numsend (2) }\
                else{\
                  if($?4 == 1) { define numsend (4) } else { define numsend 3 }
                }
		define PLANE (3)
		define WHICHLEV (0)
                #defaults
		rdnumd
                set h1=$1
                set h3='.dat'
                do ii=0,$NUMDUMPS-1,$ANIMSKIP {
                  set h2=sprintf('%04d',$ii)
                  set _fname=h1+h2+h3
                  define filename (_fname)
                  if($numsend==2){ pls  $filename $2}\
                  else{\
                   if($numsend==3){  pls  $filename $2 $3}\
                   else{\
                    if($numsend==4){ pls  $filename $2 $3 $4 $5 $6 $7}
                   }
                  }
                  #delay loop
                  #set jj=0
                  #while {jj<1} {set jj=jj+1}
		} 
		#
animpl	18	#
                if($?4 == 0) { define numsend (3) }\
                else{\
                  if($?5 == 1) { define numsend (5) } else { define numsend 4 }
                }
                #defaults
		rdnumd
                set h1=$1
                set h3='.dat'
                do ii=0,$NUMDUMPS-1,$ANIMSKIP {
                  set h2=sprintf('%04d',$ii)
                  set _fname=h1+h2+h3
                  define filename (_fname)
                  if($numsend==3){ pl  $filename $2 $3}\
                  else{\
                   if($numsend==4){  pl  $filename $2 $3 $4}\
                   else{\
                    if($numsend==5){ pl  $filename $2 $3 $4 $5 $6 $7 $8}
                   }
                  }
                  #delay loop
                  #set jj=0
                  #while {jj<1} {set jj=jj+1}
		}
animplg	18	#
                if($?4 == 0) { define numsend (3) }\
                else{\
                  if($?5 == 1) { define numsend (5) } else { define numsend 4 }
                }
                #defaults		
                set h1=$1
                do ii=0,$NUMDUMPS-1,$ANIMSKIP {
                  set h2=sprintf('%03d',$ii)
                  set _fname=h1+h2
                  define filename (_fname)
                  if($numsend==3){\
		     grdp $filename
 		     pl  0 $2 $3
		  }\
                  else{\
                   if($numsend==4){\
		            grdp $filename
		            pl  0 $2 $3 $4
		         }\
                   else{\
                    if($numsend==5){\
		                   grdp $filename
		                   pl  0 $2 $3 $4 $5 $6 $7 $8
		                }
                   }
                  }
                  #delay loop
                  #set jj=0
                  #while {jj<1} {set jj=jj+1}
		}
animplflim  15	#
                # setlimits  1.0 1.4 0 3.14159 0 1 animplflim 'dump' x1 r 0 100
                if($?4 == 0) { define numsend (3) }\
                else{\
                  if($?5 == 1) { define numsend (5) } else { define numsend 4 }
                }
                #defaults
		rdnumd
                set h1=$1
                set h3='.dat'
                do ii=0,$NUMDUMPS-1,$ANIMSKIP {
                  set h2=sprintf('%04d',$ii)
                  set _fname=h1+h2+h3
                  define filename (_fname)
                  if($numsend==3){ plflim  $filename $2 $3}\
                  else{\
                   if($numsend==4){  plflim  $filename $2 $3 $4}\
                   else{\
                    if($numsend==5){ plflim  $filename $2 $3 $4 $5}
                   }
                  }
                  #delay loop
                  #set jj=0
                  #while {jj<1} {set jj=jj+1}
		}
animplo	3	#pl	18	# pl <file> <dir> <function> <logx=1000,0000,logy=0100,0000,overlay=0010,0000,limits=0000,0001> <0 0 0 0>
                animpl $1 $2 $3 0010
animpllim 7	#
                animpl $1 $2 $3 0001 $4 $5 $6 $7
                #
animplconv 3	#
                #defaults
		# when you read in text(h1) use ''s, eg. 'dump401'
                set h1=$1
                set h3='.dat'
                do ii=15,1024,10 {
                  set h2=sprintf('%04d',$ii)
                  set _fname=h2+h1+h3
                  define filename (_fname)
		  pl $filename $2 $3
                  #delay loop
                  #set jj=0
                  #while {jj<1} {set jj=jj+1}
		}  
animpld 2	#
                #defaults
		rdnumd
		set h0='a'
                set h1='dump'
                set h3='.dat'
                do ii=0,$NUMDUMPS-1,$ANIMSKIP {
                  set h2=sprintf('%04d',$ii)
                  set _fname=h1+h2+h3
                  define filename (_fname)
                  set _fname2=h0+h1+h2+h3
                  define filename2 (_fname2)
		  #pld $filename $filename2 $1 $2
                  pld $filename adump0000.dat $1 $2
                  #delay loop
                  #set jj=0
                  #while {jj<10000} {set jj=jj+1}
		}
                # when you read in a full path for macro read do: re "fullpath"
animpldconv 3	#
                #defaults
		set h0='a'
		# h1 is dumpxxx for filename in convergence test
                set h1=$1
                set h3='.dat'
                do ii=15,1024,10 {
                  set h2=sprintf('%04d',$ii)
                  set _fname=h2+h1+h3
                  define filename (_fname)
                  set _fname2=h2+h0+h1+h3
                  define filename2 (_fname2)
		  pld $filename $filename2 $2 $3
                  #delay loop
                  #set jj=0
                  #while {jj<10000} {set jj=jj+1}
		}  
adv     1       #
		set exponent=-1
		# -1 for cartesian -2 for cyl -3 for spc
                set ii=1,$nx
                set arho=(1+.1*_time)**(exponent)*ii/ii
                set rdiff=arho-r
                set avel=.1*x/(1.0+.1*_time)
                set vdiff=avel-vx
                ctype blue 
		ptype 4 0
		points x $1
		connect x $1        
advd    1       #
		set exponent=-1
		# -1 for cartesian -2 for cyl -3 for spc
                set ii=1,$nx
                set arho=(1+.1*_time)**(exponent)*ii/ii
                set rdiff=arho-r
                set avel=.1*x/(1.0+.1*_time)
                set vdiff=avel-vx
                ctype default
		ptype 4 0
                limits x $1
                erase
                box
		points x $1
		connect x $1        
advdlim    5    #
                set ii=1,$nx
		set exponent=-1
		# -1 for cartesian -2 for cyl -3 for spc
                set arho=(1+.1*_time)**(exponent)*ii/ii
                set rdiff=arho-r
                set avel=.1*x/(1.0+.1*_time)
                set vdiff=avel-vx
                ctype default
		ptype 4 0
                limits $2 $3 $4 $5
                erase
                box
		points x $1
		connect x $1        
animadvd 1	#
		rdnumd
                set h1='dump'
                set h3='.dat'
                do ii=0,$NUMDUMPS-1,1 {
                  set h2=sprintf('%04d',$ii)
                  set _fname=h1+h2+h3
                  define filename (_fname)
                  device 0
		  plx $filename $1
                  device X11
                  advd $1
                  #delay loop
                  #set jj=0
                  #while {jj<1} {set jj=jj+1}
		}
                #
plb	2	#
                #defaults
                rdbondi
		limits $1 $2
		erase
                ctype default
		box
		ptype 4 0
		points $1 $2
		connect $1 $2
		xla $1
		yla $2
plbo	2	#
                #defaults
                rdbondi
                ctype blue
		ptype 4 0
		points $1 $2
		connect $1 $2
		xla $1
		yla $2
                  ##################
                  ##################
                  ##################
                  ##################
                  #BEGIN SUPERPOSTPROCESS MACROS
                  #
                  #
                  #
                  #
plottype  2       #
                  if(!($finaldraft)){\
                   if($1==1){\
		    ltype 0
                    ptype 4 0
                    ctype default
		   }
                   if($1==2){\
		    ltype 0
                    ptype 4 0
                    ctype red
		   }
                   if($1==3){\
		    ltype 0
                    ptype 4 0
                    ctype blue
		   }
		  }\
                  else{\
                   if($1==1){\
                    ctype default
                    ltype 0
		   }
                   if($1==2){\
                    ctype default
                    ltype 3
		   }
                   if($1==3){\
                    ctype default
                    ltype 4
		   }
		  }
invertaxes 4      #
		  set $3=1,$nx*$ny,1
		  set $4=1,$nx*$ny,1
                  do ii=0,$nx*$ny-1,1 {
 		    set $3[$nx*$ny-1-i]=$1[i]
 		    set $4[$nx*$ny-1-i]=$2[i]
  		  }
rot180del 2       # cartesion coord only
                  # used with TVDLF all vars
                  # used with ZEUS with zone centered scalars
		  set $2=1,$nx*$ny,1
                  do ii=0,$nx*$ny-1,1 {
		    set $2[$ii]=$1[$nx*$ny-1-$ii]-$1[$ii]
  		  }
                  #
rot180del2 4       # cartesion coord only
                  # used with TVDLF all vars
                  # used with ZEUS with zone centered scalars
		  set $2=1,$3*$4,1
                  do ii=0,$3*$4-1,1 {
		    set $2[$ii]=$1[$3*$4-1-$ii]-$1[$ii]
  		  }
                  #  
rot180delv 2       # cartesion coord only
                  # use with ZEUS with zone centered vectors
		  set $2=1,$nx*$ny,1
                  do ii=0,$nx*$ny-1,1 {
		    set $2[$ii]=$1[$nx*$ny-1-$ii]+$1[$ii]
  		  }
                  #
rot180dvx 2       # cartesion coord only(for LOOPF of pdump, full grid inputted)
                  #  used with ZEUS vx location as vector
                  # only thing neglected is i=-2 vx since i=N1+N1BND doesn't exist as symmetry pair
		  set $2=1,$nx*$ny,1
                  do ii=0,$nx*$ny-1,1 {
                   if( ($ii>=1)&&($ii % $nx) ){\
		    set $2[$ii]=$1[$ii]+$1[$nx*$ny-1+1-$ii]
		   }\
		   else{\
		    set $2[$ii]=0
		   }
  		  }
                  #
rot180dvx2 4       # cartesion coord only(for LOOPF of pdump, full grid inputted)
                  #  used with ZEUS vx location as vector
                  # only thing neglected is i=-2 vx since i=N1+N1BND doesn't exist as symmetry pair
		  set $2=1,$3*$4,1
                  do ii=0,$3*$4-1,1 {
                   if( ($ii>=1)&&($ii % $3) ){\
		    set $2[$ii]=$1[$ii]+$1[$3*$4-1+1-$ii]
		   }\
		   else{\
		    set $2[$ii]=0
		   }
  		  }
                  #
rot180dvy 2       # cartesion coord only(for LOOPF of pdump, full grid inputted)
                  # used with ZEUS vy location as vector
                  # only thing neglected is j=-2 vy since j=N2+N2BND doesn't exist as symmetry pair
		  set $2=1,$nx*$ny,1
                  do ii=0,$nx*$ny-1,1 {
                   if( $ii>=$nx ){\
		    set $2[$ii]=$1[$ii]+$1[$nx*$ny-1+$nx-$ii]
		   }\
		   else{\
		    set $2[$ii]=0
		   }
  		  }
                  #
rot180dvy2 4       # cartesion coord only(for LOOPF of pdump, full grid inputted)
                  # used with ZEUS vy location as vector
                  # only thing neglected is j=-2 vy since j=N2+N2BND doesn't exist as symmetry pair
		  set $2=1,$3*$4,1
                  do ii=0,$3*$4-1,1 {
                   if( $ii>=$3 ){\
		    set $2[$ii]=$1[$ii]+$1[$3*$4-1+$3-$ii]
		   }\
		   else{\
		    set $2[$ii]=0
		   }
  		  }
                  #
rot180dvc 2           # cartesion coord only(for LOOPF of pdump, full grid inputted)
                  # used with vectors on corners
                  # only thing neglected is i=-2,j=-2 v(on center) since i=N1+N1BND,j=N1+N2BND doesn't exist as symmetry pair
		  set $2=1,$nx*$ny,1
                  do ii=0,$nx*$ny-1,1 {
                   if( ($ii>=$nx+1)&&($ii % $nx) ){\
		    set $2[$ii]=$1[$ii]+$1[$nx*$ny-1-($ii-($nx+1))]
		   }\
		   else{\
		    set $2[$ii]=0
		   }
  		  }
                  #
rot180dsc 2       # cartesion coord only(for LOOPF of pdump, full grid inputted)
                  # used with scalars on corners
                  # only thing neglected is i=-2,j=-2 v(on center) since i=N1+N1BND,j=N1+N2BND doesn't exist as symmetry pair
		  set $2=1,$nx*$ny,1
                  do ii=0,$nx*$ny-1,1 {
                   if( ($ii>=$nx+1)&&($ii % $nx) ){\
		    set $2[$ii]=$1[$ii]-$1[$nx*$ny-1-($ii-($nx+1))]
		   }\
		   else{\
		    set $2[$ii]=0
		   }
  		  }
                  #
rot180diff1 2     # cartesion coord only (for LOOPF of pdump, full grid inputted)
                  # NOT CURRENTLY USED
		  set $2=1,$nx*$ny,1
 		  set ii=0,$nx*$ny-1,1
		  set tempvar=$1 if( (ii>=$nx+1)&&(ii % $nx) )
                  set jj=0      
                  do ii=0,$nx*$ny-1,1 {
                   if( ($ii>=$nx+1)&&($ii % $nx) ){\
		    set $2[$ii]=tempvar[jj]+tempvar[($nx-1)*($ny-1)-1-jj]
                    set jj=jj+1
		   }\
		   else{\
		    set $2[$ii]=0
		   }
  		  }
                  #
writedump  1         # can use this to output any vectors/scalars and interp them using pp run
		     define print_noheader (1)
		     #read {r$!!ext 1 en$!!ext 2 pot$!!ext 3 v$!!extx 4 v$!!exty 5 v$!!extz 6 b$!!extx 7 b$!!exty 8 b$!!extz 9}
                     set it=sprintf('%10s\n','\# fileversion filetype')
                     set crap1=1
                     set crap2=$DTYPE
		     print $1 '%6s%6d %6d\n' {it crap1 crap2}
                     set it=sprintf('%6s\n','\# t SAMPLE ZONEC')
                     set crap1=$time
                     set crap2=0
                     set crap3=0
		     print + $1 '%6s%6d %6d %6d\n' {it crap1 crap2 crap3}
                     set it3=sprintf('%6s\n','\# nx ny nz')
		     set crap1=$nx
		     set crap2=$ny
		     set crap3=1
		     print + $1 '%6s%6d %6d %6d\n' {it3 crap1 crap2 crap3}
                     set it=sprintf('%s\n','\# rho u pot  vx1  vx2 vx3 bx1 bx2 bx3')
		     print + $1 '%s' {it}
                     print + $1 '%21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g\n' {r2d en2d pot2d v2dx v2dy v2dz b2dx b2dy b2dz}
                     #print + $1 '%21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g\n' {r en pot vx vy vz bx by bz}
		     define print_noheader (0)
                     #
thetaphiavg   5      # (e.g. thetaphiavg PI/6 infun outfun outx1 0)
		    # average is over all phi
                     set $3=1,$nx,1
		     set $3=0*$3
		     set $4=1,$nx,1
		     set $4=0*$4
		     set rangex2x3=1,$nx,1
		     set rangex2x3=0*rangex2x3
		     do ii=0,$nx*$ny*$nz-1,1{
		      set indexi=INT($ii%$nx)
		      set indexj=INT(($ii%($nx*$ny))/$nx)
		      set indexk=INT($ii/($nx*$ny))
		      if( ABS(x2[$ii]-PI/2)<$1 ){\
		       # full theta-phi differential term
		       set $3[indexi]=$3[indexi]+$2[$ii]*sin(x22[$ii])*dx22[$ii]*dx32[$ii]
		       # just value
		       #set $3[indexi]=$3[indexi]+$2[$ii]
		       set $4[indexi]=x12[$ii]
		       set rangex2x3[indexi]=rangex2x3[indexi]+sin(x22[$ii])*dx22[$ii]*dx32[$ii]
		      }
		     }
		     if($5==0){\
		      #if(2*$1>(x2[$nx*$ny-1]-x2[0])){ set rangex2=(x2[$nx*$ny-1]-x2[0]) } else{ set rangex2=2*$1 }
		      #set rangex3=x3[$nx*$ny*$nz-1]-x3[0]
		      # correct for degeneracies
		      #if($ny==1){ set rangex2=dx22[0] }
		      #if($nz==1){ set rangex3=dx32[0] }		     
 		      #set $3=$3/(rangex2*rangex3)
 		      set $3=$3/(rangex2x3)
		     }
		     #
accelerations  2     #
		     # accelerations ($d\mathbf{v}/dt$)
		     #
		     # define whether using limiter or not
		     #
		     # lagrangian derivative term: (v\cdot\nalbda)v
		     echo "vdotgradv"
		     vdotgradv v v atv $1 $2
		     set atvx=-atvx
		     set atvy=-atvy
		     set atvz=-atvz
		     #
		     # gravitational acceleration
		     echo "grad pot"
		     grad pot aphi $1 $2
		     set aphix=-aphix
		     set aphiy=-aphiy
		     set aphiz=0
		     #
		     # pressure gradient acceleration
		     set p=($gam-1)*en
		     echo "grad p"
		     grad p ap $1 $2
		     set apx=-apx/r
		     set apy=-apy/r
		     set apz=0
		     #
		     # magnetic pressure acceleration
		     #
		     echo "grad b2"
		     set b2=bx**2+by**2+bz**2
		     grad b2 abp $1 $2
		     set abpx=-0.5*abpx/r
		     set abpy=-0.5*abpy/r
		     set abpz=0
		     #
		     # magnetic tension acceleration
		     echo "grad bdotgradb"
		     vdotgradv b b atb $1 $2
		     set atbx=atbx/r
		     set atby=atby/r
		     set atbz=atbz/r
                     #
                     # the above has 3 components!
		     set amagx=atbx+abpx
		     set amagy=atby+abpy
		     set amagz=atbz+abpz
		     # now look at ratios
		     set rtvx=atvx/aphix
		     set rtvy=atvy/aphix
		     set rtvz=atvz/aphix		     
		     set rpx=apx/aphix
		     set rpy=apy/aphix		     
		     set rbpx=abpx/aphix
		     set rbpy=abpy/aphix
		     set rtbx=atbx/aphix
		     set rtby=atby/aphix
		     set rtbz=atbz/aphix
		     set totaccx=atvx+aphix+apx+abpx+atbx
		     set totaccy=atvy+aphiy+apy+abpy+atby
		     set totaccz=atvz+aphiz+apz+abpz+atbz
		     set totacc=sqrt(totaccx**2+totaccy**2+totaccz**2)
		     set forcex=totaccx*r
                     set forcey=totaccy*r
                     set forcez=totaccz*r
		     set rtotx=totaccx/aphix
		     set rtoty=totaccy/aphix
		     set rtotz=totaccz/aphix
		     #
b2time	1	#
		rdnumd
		set b2vstime=1,$NUMDUMPS,1
		set fieldtime=1,$NUMDUMPS,1
                set h1=$1
                set h3='.dat'
                do ii=0,$NUMDUMPS-1,$ANIMSKIP {
                  set h2=sprintf('%04d',$ii)
                  set _fname=h1+h2+h3
                  define filename (_fname)
		  rd $filename
		  polefield 1
		  set b2vstime[$ii]=avgmax
		  set fieldtime[$ii]=$time
		}
revertmin1 0    #
		rd 0_loss.dat
		set ii=0,dimen(min1)-1,1
		set min1new=min1 if(!(ii%(INT(_DTd/_DTloss))))
		der fieldtime min1new td min1d
plotcontour 0    #
                 da data.txt
		 read {x 1 y 2 z 3}
                 # size of data
                 define rnx 6
                 define rny 6
                 # limits of data
                 define txl (0.0)
                 define txh (100.0)
                 define tyl (0.0)
                 define tyh (100.0)
		 # uniform grid dx for each direction
		 set temp=($txh-$txl)/($rnx-1)
                 define realdx (temp)
		 set temp=($tyh-$tyl)/($rny-1)
                 define realdy (temp)
                 # number of contours
                 define cres (10)
                 # whether to have solid contours or shaded for negative values
                 define SOLIDCONTOURS (0)
                 #
		 set ii=1,$rnx*$rny,1
		 set x=(ii - $rnx*int((ii-0.5)/$rnx) - 0.5)*$realdx+$txl
		 set y=((int((ii-0.5)/$rnx) + 0.5)*$realdy)+$tyl
		 set ix = int((x -$txl)/$realdx)
		 set iy = int((y -$tyl)/$realdy)
                 limits $txl $txh $tyl $tyh
 		 image ($rnx,$rny) $txl $txh $tyl $tyh
                 set image[ix,iy] = z
                 erase
		 minmax min max
                 echo "min:"$min "max:"$max
		 	if(($SOLIDCONTOURS==0)&&($min*$max < 0.)) {\
                        echo "pt:1"
			define delta (($max-$min)/$cres)
                        echo "delta:"$delta
                        if($min>-$delta){\
                          define min -$delta
                        }
			set lev=$min,-$delta,$delta
			levels lev
			ltype 4
			contour
			#
                        if($max<$delta){\
                          define max $delta
                        }
			set lev=$delta,$max,$delta
			levels lev
			ltype 0
			contour
		}
		if(($min<=0)&&($max<=0)) {\
                        echo "pt:2"
			set lev=$min,$max,($max-$min)/$cres
			levels lev
			#ltype 2
                        ltype 0
			contour
		}
		if(($SOLIDCONTOURS==1)||($min>=0)&&($max>=0)) {\
                        echo "pt:3"
			set lev=$min,$max,($max-$min)/$cres
			levels lev
                        ltype 0
			contour
                }
                box
                #
 #      set radius=sqrt(x12**2+x22**2+x32**2)+1E-31
 # 	set ndot=4*pi*radius**2 * r *(vx*x12/radius+vy*x22/radius+vz*x32/radius)
 # 	set dmx=r*vx*dx2*dx3
 # 	set dmy=r*vy*dx1*dx3
 # 	set dmx=r*vz*dx1*dx2
 #		read {m0 1 m1 2 m2 3 m3 4 m4 5 m5 6 m6 7 m7 8 m8 9 m9 10 m10 11 m11 12 m12 13 m13 14 m14 15 \
 #		m15 16 m16 17 m17 18 m18 19 m19 20 m20 21 m21 22 m22 23 m23 24 m24 25 m25 26 m26 27 m27 28 m28 29}
 # ((gam-1.0)*s[2][kk][jj][ii]+s[1][kk][jj][ii]*s[3][kk][jj][ii])*g[2][2][ii]*g[2][3][ii]*g[2][4][jj]*dx[1][2][jj]*dx[1][3][kk];
	#
animgpls	17	# animgpls 'dump' rho
                if($?3 == 0) { define numsend (2) }\
                else{\
                  if($?4 == 1) { define numsend (4) } else { define numsend 3 }
                }
                #defaults
		#
		define NUMDUMPS (11)
                set h3='.0000'
                set h1=$1
                do ii=0,$NUMDUMPS-1,$ANIMSKIP {
                  set h2=sprintf('%03d',$ii)
                  set _fname=h1+h2+h3
                  define filename (_fname)
		  #re twod.m
		  rdp $filename
		  #smstart
		  define filename (0)
                  if($numsend==2){ pls  $filename $2}\
                  else{\
                   if($numsend==3){  pls  $filename $2 $3}\
                   else{\
                    if($numsend==4){ pls  $filename $2 $3 $4 $5 $6 $7}
                   }
                  }
                  #delay loop
                  #set jj=0
                  #while {jj<1} {set jj=jj+1}
		} 
                #
grey_sincos     # draw a grey scale image of a sin(x)cos(y) surface
		erase
		define nx (128)
		define ny (128)
                image ( $nx , $ny ) 0 1 0 1
                glevels do(-1,1,0.05)
                set x local set y local
                set y=0,$ny-1
                do x=0,$nx{
                   set image[$x,y] = sin(2*$x*2*pi/$nx)*cos(5*y*2*pi/$ny)
                }
                lim 0 1 0 1
                box greyscale 40 40 8
                define 1 1
                define 1 ? {Now draw that using dither; hit <CR> when ready}
                ERASE
                dither x y -1.1 2 5
                set x=x+random(dimen(x))*2e-3
                set y=y+random(dimen(y))*2e-3
                define ptype local define ptype |
                ptype 1 1
                box poi x y
                ptype $ptype
                #
maskit    2    #
		set mask=0,$nx*$ny*$nz-1,1
		do ii=0,$nx*$ny*$nz-1,1{
		   set indexi=INT($ii%$nx)
		   set indexj=INT(($ii%($nx*$ny))/$nx)
		   set indexk=INT($ii/($nx*$ny))
		   
		   if(abs(indexi-$nx/2)<$2){ set insidex=1} else{ set insidex=0 }
		   if(abs(indexj-$ny/2)<$2){ set insidey=1} else{ set insidey=0 }
		   if(abs(indexk-$nz/2)<$2){ set insidez=1} else{ set insidez=0 }
		   # want the x,y,z condition to not fucking matter if nx,ny,nz=1
		   if($nx==1) { set nxis1=1 } else { set nxis1=0 }
		   if($ny==1) { set nyis1=1 } else { set nyis1=0 }
		   if($nz==1) { set nzis1=1 } else { set nzis1=0 }
		   #
		   if((nxis1||insidex)&&(nyis1||insidey)&&(nzis1||insidez)){ set $1[$ii]=0 } else {set $1[$ii]=1 }
		   
                }
parperpparts 4   # parperparts a b ab 2
		set $3x=$2x*0
                set $3y=$2y*0
                set $3z=$2z*0
		if($4==2){ set $3x=($1x*$2x+$1y*$2y)/sqrt($2x*$2x+$2y*$2y) }
		if($4==2){ set $3y=($1x*$2y-$1y*$2x)/sqrt($2x*$2x+$2y*$2y) }
		if($4==3){ set $3x=($1x*$2x+$1y*$2y+$1z*$2z)/sqrt($2x*$2x+$2y*$2y+$2z*$2z) }
		# just consistency check
		if($4==4){ set $3x=($1x*$2x+$1y*$2y+$1z*$2z)/(sqrt($2x*$2x+$2y*$2y+$2z*$2z)*sqrt($1x*$1x+$1y*$1y+$1z*$1z)) }
		if($4==3){ set $3y=0 }
		if($4==3){ set $3z=0 }
		#
dtcalc   0  #		
		        set cs2=$ccs/(r*$ccs/($gam*($gam-1)*en)+1)
		        set va2=$cca/(r*$cca/(bx**2+by**2+bz**2)+1)
		        set dtinvx=(vx+sqrt(va2+cs2))/dx1
		        set dtinvy=(vy+sqrt(va2+cs2))/(g2*dx2)		
crapit	0	#
	#
