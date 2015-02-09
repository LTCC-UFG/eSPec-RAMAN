#!/bin/bash
#  Raman script for eSPec
#
#  This is a simple script to run RIXS calculations with eSPec
#
#   
#  required programs:
#  eSPec - https://github.com/LTCC-UFG/eSPec
#  raman - https://github.com/LTCC-UFG/eSPec-RAMAN
#
#
# Vinicius Vaz da Cruz - viniciusvcruz@gmail.com
#
# Goiania, 21st of January of 2015
#

#eSPec path
espec=/home/vinicius/eSPec_v0.7/espec_v07.x
#raman-eSPec path
raman=/home/vinicius/programming/LTCC/eSPec-RAMAN/raman

#----------- Script modes ---------------#
runtype=$1

#
# -all   runs all three steps
# 
# -init  runs only the initial propagation. i.e. propagation of |0> on the core-excited potential
# 
# -cond  olnly generates the |Phi(0)> intermediate wavepackets (-init must have been run previously)
# 
# -fin   runs only the final propagation. i.e. propagation of |Phi(0)> on the final potential 
#        (-init must have been run previously)
# 
# -cfin  combination of -cond and -fin
# 

#-----------General input parameters---------------#
input=$2
# short name to be the base file name for input, output and result files
jobid=`grep -i -w jobid $input | awk '{printf $2}'`
# dimension .1D .2D .2DCT
dim=`grep -i -w dimension $input | awk '{printf $2}'`
# number of discretization points used. for .2D use npoints="n1 n2"
work=`grep -i -w npoints $input | awk '{printf $1}'`
npoints=`grep -i -w npoints $input | sed "s/\<$work\>//g"`
#-------------------------
#
initial_wf=`grep -i -w initial_wf $input | awk '{printf $2}'`
if [ "$initial_wf" -eq ".CALC" ] || [ -z "$initial_wf" ]; then
    initial_pot=`grep -i initial_pot $input | awk '{printf $2}'`
    mode='.CALC'
else
    initial_pot=`grep -i -w initial_wf $input | awk '{printf $2}'`
    mode='.GETC'
fi
#
# potential files---------
decaying_pot=`grep -i decaying_pot $input | awk '{printf $2}'`
final_pot=`grep -i final_pot $input | awk '{printf $2}'`
#mass in amu for .2D just use mass="value1 value2" .2DCT mass="value1 value2 value3"
work=`grep -i -w mass $input | awk '{printf $1}'`
mass=`grep -i -w mass $input | sed "s/\<$work\>//g"`
#core excited lifetime
Gamma=`grep -i gamma $input | awk '{printf $2}'`
# time step used in the propagation
step=`grep -i step $input | awk '{printf $2}'`
#values of the detuning desired
work=`grep -i detun $input | awk '{printf $1}'`
all_detunings=`grep -i detun $input | sed "s/\<$work\>//g"`
#all_detunings="-2.0 -1.0 -0.1 0.0 0.1 1.0 2.0 4.0"

#recomended propagation time based on Gamma
decay_threshold='1e-3'
init_time=$(awk "BEGIN {print (-log($decay_threshold)/($Gamma/27.2114))/41.3411}")

# propagation time on final state
fin_time=`grep -i fin_time $input | awk '{printf $2}'`

#---------------Initial Propagation---------------#

echo "-----------------"
echo "eSPec-Raman script"
echo "------------------"

echo
echo "job running on" `hostname`
date
echo

ulimit -s unlimited


if [ "$runtype" == "-all" ] || [ "$runtype" == "-init" ]; then

    echo ' Starting initial propagation'
    echo
    echo "propagation time on decaying potential: $init_time fs"


    cat $initial_pot > pot.inp
    cat $decaying_pot >> pot.inp

    cat > input.spc <<EOF
*** eSPec input file ***
========================
**MAIN
*TITLE
 +++++ $jobid Raman init +++++
*DIMENSION
$dim
 $npoints/
*POTENTIAL
.FILE
pot.inp
*MASS
$mass/
*TPCALC
.PROPAGATION
.TD
*INIEIGVC
$mode
*CHANGE
.YES
*PRTCRL
.PARTIAL
*PRTPOT
.YES
*PRTEIGVC
.YES
*PRTVEFF
.NO
*PRTEIGVC2
.YES


**TI
*TPDIAG
.MTRXDIAG 
*NIST
 5  0/
*ABSTOL
 1D-6/

**TD
*PROPAG
.PPSOD
 0.0  $init_time  $step/
*PRPGSTATE
 0/
*TPTRANS
.ONE
*PRPTOL
 3.0D0/
*NPROJECTIONS
 6  1
*FOURIER
 14/

**END
EOF

    time $espec > ${jobid}_init.out
    mkdir wf_data
    mv eigvc_*.dat ReIm_*.dat wf_data/

    echo 'Initial propagation done!'
    echo
    #-------------------------------------------------#

fi


if [ "$runtype" == "-all" ] || [ "$runtype" == "-cond" ] || [ "$runtype" == "-cfin" ]; then
    #---------------|Phi(0)> calculation--------------#
    echo 'Generating initial conditions for second propagation'
    echo


    if [ "$dim" == ".1D" ]; then

	Vg_min=`sed '/#/ d' $initial_pot | awk '{printf $2"\n"}' | perl -e 'print sort { $a<=>$b } <>' | head -1`
	#echo "ground state min" $Vg_min
	echo "ground state min" $Vg_min > $jobid.log
	pos=`cat -n $initial_pot | grep $Vg_min | awk '{printf $1"\n"}' | head -1`
	#echo "ground state minimum position" $pos
	echo "ground state minimum position" $pos >> $jobid.log
	Vd_min=`sed '/#/ d' $decaying_pot | awk '{printf $2"\n"}' | perl -e 'print sort { $a<=>$b } <>' | head -1`
	#echo "decaying state min" $Vd_min
	echo "decaying state min" $Vd_min >> $jobid.log
	Vd=`cat -n $decaying_pot | grep -w " $pos" | awk '{printf $3}'`
	#echo "decaying state vertical" $Vd
	echo "decaying state vertical" $Vd >> $jobid.log
	Vf_min=`sed '/#/ d' $final_pot | awk '{printf $2"\n"}' | perl -e 'print sort { $a<=>$b } <>' | head -1`
	#echo "final state min" $Vf_min
	echo "final state min" $Vf_min >> $jobid.log

	if [ -f "${jobid}_init.out" ]; then

	    E0=`grep '|     0        |' ${jobid}_init.out | awk '{printf $4}'`
	    echo "Initial energy" $E0
	    echo "Initial energy" $E0 >> $jobid.log

	else

	    echo "failed to find initial propagation output file ${jobid}_init.out"
	    echo "Attempting to read E0 from input file"
	    E0=`grep -i -w "E0" | awk '{printf $2}'`
	    if [ -z "$E0" ]; then
		echo "could not find the value of E0"
		echo "please check your input"
		exit 666
	    fi

	fi

	omres=$(awk "BEGIN {print $Vd - $Vg_min - $E0}")
	echo "resonance frequency: $omres"
	echo "resonance frequency: $omres" >> $jobid.log
	Eres=$(awk "BEGIN {print $omres - $Vd_min + $E0}")
	echo "shifted resonance frequency: $Eres"
	echo "shifted reson. frequency: $Eres" >> $jobid.log
	Vgf_gap=$(awk "BEGIN {print $Vf_min - Vg_min}")
	echo "ground to final gap: $Vgf_gap"
	echo "ground to final gap: $Vgf_gap" >> $jobid.log

    elif [ "$dim" == ".2D" ]; then

	Vg_min=`grep -i Vg_min $input | awk '{printf $2}'`
	Vd_min=`grep -i Vd_min $input | awk '{printf $2}'`
	Vf_min=`grep -i Vf_min $input | awk '{printf $2}'`
	Vd=`grep -i Vd_vert $input | awk '{printf $2}'`

	E0=`grep '|     0        |' ${jobid}_init.out | awk '{printf $4}'`
	echo "Initial energy" $E0
	echo "Initial energy" $E0 >> $jobid.log
	omres=$(awk "BEGIN {print $Vd - $Vg_min - $E0}")
	echo "resonance frequency: $omres"
	echo "resonance frequency: $omres" >> $jobid.log
	Eres=$(awk "BEGIN {print $omres - $Vd_min + $E0}")
	echo "shifted resonance frequency: $Eres"
	echo "shifted reson. frequency: $Eres" >> $jobid.log
	Vgf_gap=$(awk "BEGIN {print $Vf_min - Vg_min}")
	echo "ground to final gap: $Vgf_gap"
	echo "ground to final gap: $Vgf_gap" >> $jobid.log

    fi

    nfiles=`ls wf_data/ReIm_*.dat | awk '{printf $1"\n"}' | tail -1 | cut -c 14-17`
    last_file=`ls wf_data/ReIm_*.dat | awk '{printf $1"\n"}' | tail -1`
    rtime=`head -1 $last_file | awk '{printf $3}'`

    echo
    echo "number of wavepacket files: $nfiles"
    #ndetun=`echo $all_detunings | wc -w`

    cat > raman.inp <<EOF
# eSPec-Raman input

*Main
dimension
$dim
$npoints
mass: $mass
filename: wf_data/ReIm_
nfiles: 1 $nfiles
timeinterval: 0.000 $rtime

*Propagation
Ereso: $Eres
detuning: $all_detunings
Fourier: 10
Window
.EXPDEC $Gamma
EOF

    time $raman > ${jobid}_raman.out

    echo
    echo 'initial conditions generated!'
    echo

    #-------------------------------------------------#
fi

if [ "$runtype" == "-all" ] || [ "$runtype" == "-fin" ] || [ "$runtype" == "-cfin" ]; then

    #---------------Final Propagation-----------------#
    echo ' Starting final propagation'

    for detun in `echo $all_detunings`
    do
	
	echo "running detuning = $detun"

	detun2=${detun}00000
	file=wp_$detun2.dat

	cat $final_pot >> $file

	cat > input.spc <<EOF
*** eSPec input file ***
========================
**MAIN
*TITLE
 +++++ $jobid Raman fin +++++
*DIMENSION
$dim
 $npoints/
*POTENTIAL
.FILE
$file
*MASS
$mass/
*TPCALC
.SPECTRUM
.TD
*INIEIGVC
.GETC
*CHANGE
.YES
*PRTCRL
.PARTIAL
*PRTPOT
.NO
*PRTEIGVC
.NO
*PRTVEFF
.NO
*PRTEIGVC2
.NO

**TD
*PROPAG
.PPSOD
 0.0  $fin_time  $step/
*PRPGSTATE
 0/
*TPTRANS
.ONE
*PRPTOL
 3.0D0/
*NPROJECTIONS
 6  1

**SPECTRUM
*FOURIER
 14/
*WINDOWING
.SG
1D-5/

**END
EOF

	time $espec > ${jobid}_$detun.out
	sed -n "/Spectrum:/,/Spectrum done./p" ${jobid}_$detun.out | sed "/Spectrum/ d" | sed "/*/ d" | sed "/=/ d" > temp
	cat temp | awk '{printf $1" "$2"\n"}' > bkp-${jobid}_$detun.spec
	rm temp

	echo "propagation for detuning = $detun done!"
	echo


	num=`cat -n  bkp-${jobid}_$detun.spec | tail -1 | awk '{printf $1}'`
	norm=`head -1 $file | awk '{printf $8}'`
	omres=`grep "resonance frequency:" $jobid.log | awk '{printf $3}'`
	omres=$(awk "BEGIN {print $omres * 27.2114}")
	Vgf_gap=`grep "ground to final gap:" $jobid.log | awk '{printf $5}'`
	E0=`grep "Initial energy" $jobid.log | awk '{printf $3}'`
	Vgf_gap=$(awk "BEGIN {print ($Vgf_gap - $E0) * 27.2114}")
	echo "$num $omres $detun $Vgf_gap $norm bkp-${jobid}_$detun.spec " > sinp
	./sshift < sinp > ${jobid}_$detun.spec

    done

    #-------------------------------------------------#

fi

echo
echo 'eSPec-Raman script finished'
