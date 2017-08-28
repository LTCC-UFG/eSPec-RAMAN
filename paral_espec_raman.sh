#!/bin/bash
#  Raman script for eSPec
#
#  This is a script to run RIXS calculations with eSPec
#  for the specific case of our 2D + 1D model
#   
#
#  required programs:
#  eSPec - https://github.com/LTCC-UFG/eSPec
#  raman - https://github.com/LTCC-UFG/eSPec-RAMAN
#
#
# Vinicius Vaz da Cruz - viniciusvcruz@gmail.com
#
# Goiania, 27th of January of 2015
#

#%header1%
#%header2%
#%header3%

#eSPec path
espec=%especpath%
#raman-eSPec path
raman=%ramanpath%
#fcorrel path
fcorrel=%fcorrelpath%

#----------- Script modes ---------------#
runtype=$1

#
# -all   runs all three steps
# 
# -init  runs only the initial propagation. i.e. propagation of |0> on the core-excited potential
#
# -fc    only computes the bending franck-condon factors
# 
# -cond  only generates the |Phi(0)> intermediate wavepackets (-init must have been run previously)
# 
# -fin   runs only the final propagation. i.e. propagation of |Phi(0)> on the final potential 
#        (-init must have been run previously)
# 
# -cfin  combination of -cond and -fin
#
# -cross computes only the final cross section (all other steps must have been run previously)
#
# -clean deletes all files generated during the calculation
# 
#
############################3
#
# -xas X-ray Absorpotion Spectrum considering the 2D+1D model (this option also runs the -init and -fc steps from above)
#
# -xascs computes only the final part to get the XAS spectrum (you must have run -all, or -init and -fc previously)
#
####################################
#
# -self REXS cross section with self-absorption factor ( you must have manually run -all and -xas previously!!)
#



#---------------
if [ "$runtype" == "-clean" ]; then


    rm intens*.dat inwf*.dat wp-vc*.dat wp-vc*.inp
    rm -r fin_vc* *cond-vc* wf_data/
    rm -r fc_*.dat fc*.out fc*.inp fcond.dat
    rm -r debug*.dat veff_0001.dat movie.gplt pot.inp bpot.inp
    rm -r fcorrel_vc*.dat fcorrel*.dat correl.inp
    rm *.log *~ *.out
    rm xas-fcorrel.dat
    rm input.spc initial_spc.aux
    rm *.spec raman.inp
    exit
fi

echo "-----------------"
echo "eSPec-Raman script"
echo "------------------"

echo
echo "job running on" `hostname`
date
echo

ulimit -s unlimited

#-----------General input parameters---------------#
input=$2
# model used in the calculation
model=`grep -i model $input | awk '{printf $2}'`
if [ -z "$model" ] || [ "$model" == "2d+1d" ] || [ "$model" == "1d+1d" ] || [ "$model" == "nd+1d" ]; then
    tpmodel=0
elif [ "$model" == "2d" ] || [ "$model" == "1d" ] || [ "$model" == "nd" ] || [ "$model" == "simple" ] || [ "$model" == "pure" ]; then
    tpmodel=1
elif [ "$model" == "1d+2d" ] || [ "$model" == "1d+3d" ] || [ "$model" == "1d+md" ] || [ "$model" == "2d+md" ] || [ "$model" == "nd+md" ]; then
    tpmodel=2
else
    echo "error: Invalid model choice $model"
    exit
fi

echo
if [ $tpmodel -eq 0 ]; then
    echo 'This calculations uses the nD(TD) + 1D(TI) model'
elif [ $tpmodel -eq 1 ]; then
    echo 'This calculation uses the simple nD(TD) model'
elif [ $tpmodel -eq 2 ]; then
    echo 'This calculations uses the nD(TD) + mD(TI) model'
fi
echo

# short name to be the base file name for input, output and result files
jobid=`grep -i jobid $input | awk '{printf $2}'`
# dimension .1D .2D .2DCT
dim=`grep -i -w dimension $input | awk '{printf $2}'`
# if .2DCT is inputed, then we need the value of the cross term cos(\theta)
if [ "$dim" == ".2DCT" ]; then
    cross=`grep -i -w cross_term $input | awk '{printf $2}'`
    crosskey="*CROSS_TERM"
    if [ -z "$cross" ]; then
	echo "dimension .2DCT, but failed to read cos(\theta) cross term value"
	echo "please, check your input"
	exit 666
    fi
fi

# number of independent modes to be included in the FC factor calculations
dim_fc=`grep -i fc_dimens $input | awk '{printf $2}'`
if [ -z "$dim_fc" ]; then
    echo "default single mode FC factors."
    dim_fc='.1D'
else
    echo "$dim_fc independent modes will be included in the FC factors."
fi

# number of discretization points used. for .2D use npoints="n1 n2"
work=`grep -i -w npoints $input | awk '{printf $1}'`
npoints=`grep -i -w npoints $input | sed "s/\<$work\>//g"`
#-------------------------
#
initial_wf=`grep -i -w initial_wf $input | awk '{printf $2}'`
if [ "$initial_wf" == ".CALC" ] || [ -z "$initial_wf" ]; then
    initial_pot=`grep -i initial_pot $input | awk '{printf $2}'`
    mode='.CALC'

    chkst=`grep -i -w init_state $input | awk '{printf $1}'`
    if [ -z "$chkst" ]; then
	instate='0'
	nist='10'
    else
	instate=`grep -i -w init_state $input | awk '{printf $2}'`
	if [ $instate -gt 10 ]; then
	    nist=$(echo 2 + $instate | bc)
	    echo "number of vibrational eigenstates to be calculated changed to $nist"
	fi
    fi

else
    initial_pot=`grep -i -w initial_wf $input | awk '{printf $2}'`
    mode='.GETR'
    instate='0'
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
work=`grep -i -w detuning $input | awk '{printf $1}'`
all_detunings=`grep -i -w detuning $input | sed "s/\<$work\>//g"`
#or alternatively values of omega desired
work=`grep -i -w omega_incoming $input | awk '{printf $1}'`
all_omega=`grep -i -w omega_incoming $input | sed "s/\<$work\>//g"`
#all_detunings="-2.0 -1.0 -0.1 0.0 0.1 1.0 2.0 4.0"

#recomended propagation time based on Gamma
decay_threshold='1e-3'
recom_init_time=$(awk "BEGIN {print (-log($decay_threshold)/($Gamma/27.2114))/41.3411}")
init_time=`grep -i -w init_time  $input | awk '{printf $2}'`

# propagation time on final state
fin_time=`grep -i fin_time $input | awk '{printf $2}'`

# absorbing conditions
absorb_cond=`grep -i -w absorb_cond $input | awk '{printf $1}'`
#bug this was inverted
if [ -z "$absorb_cond" ]; then 
    abs=" "
    abs1=" "
    abstren=" "
    absrang=" "
else
    abs="*ABSORBING"
    abs1=".SMOOTHW"
    work=`grep -i -w absorb_cond $input | awk '{printf $1}'`
    abstren=`grep -i -w absorb_cond $input | sed "s/\<$work\>//g"`"/"
    work=`grep -i -w absorb_range $input | awk '{printf $1}'`
    absrang=`grep -i -w absorb_range $input | sed "s/\<$work\>//g"`"/"
fi

# fft supergaussian window of the final spectrum
window=`grep -i -w window $input | awk '{printf $1}'`
if [ -z "$window" ]; then
    window="1e-5"
else
    window=`grep -i -w window $input | awk '{printf $2}'`
fi


#shift_spectrum
doshift=`grep -i -w shift $input | awk '{printf $1}'`
if [ -z "$doshift" ]; then
    doshift="n"
else
    doshift="y"
fi

#print level
print_level=`grep -i -w print_level $input | awk '{printf $2}'`
if [ -z "$print_level" ]; then
    print_level="essential"
fi

#transition dipole moments
dipole_id=`grep -i -w dipole_id $input | awk '{printf $2}'`
if [ -z "$dipole_id" ]; then
    dip_id="n"
else
    dip_id="y"
    echo "dipoles from initial state to decaying state will be ignored (not implemented)"
fi

dipole_df=`grep -i -w dipole_df $input | awk '{printf $2}'`
if [ -z "$dipole_df" ]; then
    dip_df="n"
else
    dip_df="y"
fi

#---------------Initial Propagation---------------#


if [ "$runtype" == "-all" ] || [ "$runtype" == "-init" ] || [ "$runtype" == "-xas" ]; then

    echo 
    echo ' Starting initial propagation'
    echo

    checktime=`echo "$init_time < $recom_init_time" | bc -l`
    if [ -z "$init_time" ]; then
	echo "recomended initial propagation time will be used!"
	init_time=`echo $recom_init_time`
    elif [ "$checktime" -eq "1" ]; then
	echo "WARNING: initial propagation time provided is smaller than the recomended one!"
	echo "changing initial propagation time from $init_time to $recom_init_time!"
	init_time=`echo $recom_init_time`
    else
	echo "user provided initial propagation time will be used!"
    fi

    echo "propagation time on decaying potential: $init_time fs"

    echo "vibrational state to be propagated : $instate"

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
$crosskey
$cross/

**TI
*TPDIAG
.MTRXDIAG 
*NIST
 12  0/
*ABSTOL
 1D-6/

**TD
*PROPAG
.PPSOD
 0.0  $init_time  $step/
*PRPGSTATE
 $instate/
*TPTRANS
.ONE
*PRPTOL
 3.0D0/
*NPROJECTIONS
 6  1
*FOURIER
 14/
$abs
$abs1
$abstren
$absrang

**END
EOF
    time $espec > ${jobid}_init.out

    #cp ReIm_0001.dat  st_0.dat

    if [ -d wf_data ]; then
	echo "previous wf_data will be replaced"
    else
	mkdir wf_data
    fi
    mv eigvc_*.dat ReIm_*.dat wf_data/

    echo 'Initial propagation done!'
    echo
    #-------------------------------------------------#
    
    
    echo 'saving XAS correlation function to file xas-fcorrel.dat'
    if [ "$dim" == ".1D" ]; then
	sed -n "/Auto-correlation function/,/ Eigenvector propagated/p"  ${jobid}_init.out | sed "/Eig/ d" | sed "/==/ d" | sed "/t\/fs/ d" | sed "/Auto/ d" | sed '/^\s*$/d' | awk '{printf $1" "$4" "$5"\n"}' > xas-fcorrel.dat
    else
	sed -n "/Partial auto-correlation function/,/ Eigenvector propagated/p"  ${jobid}_init.out | sed "/Eig/ d" | sed "/==/ d" | sed "/t\/fs/ d" | sed "/Partial/ d" | sed '/^\s*$/d' | awk '{printf $1" "$4" "$5"\n"}' > xas-fcorrel.dat
    fi
    #--------------------------------------------------#

    if [ "$print_level" == "minimal" ]; then
	rm input.spc pot.inp
	rm wf_data/eigvc_* movie.gplt veff_0001.dat
    elif [ "$print_level" == "essential" ]; then
	rm input.spc pot.inp 
	rm wf_data/eigvc_* movie.gplt veff_0001.dat
    elif [ "$print_level" == "intermediate" ]; then
	rm input.spc pot.inp movie.gplt veff_0001.dat 
    elif [ "$print_level" == "savedisk" ]; then
	rm input.spc pot.inp
	rm wf_data/eigvc_*
    fi

fi

#---------------SPECIAL CASE: 1D Franck-Condon--------------#
if [ "$runtype" == "-all" ] || [ "$runtype" == "-fc" ] || [ "$runtype" == "-xas" ] && [ $tpmodel -eq 0 ] && [ "$dim_fc" == ".1D" ]; then

    # potential files
    bendinit_pot=`grep -i -w bend_init_pot $input | awk '{printf $2}'`
    benddecay_pot=`grep -i -w bend_decay_pot $input | awk '{printf $2}'`
    bendfin_pot=`grep -i -w bend_fin_pot $input | awk '{printf $2}'`

    bend_np=`grep -i -w bend_npoints $input | awk '{printf $2}'`
    bend_mass=`grep -i -w bend_mass $input | awk '{printf $2}'`

    echo 'Computing Franck-Condon Factors'
    echo

    # bending states variables
    nvc=`grep -i -w nvc $input | awk '{printf $2}'`
    echo "number of vc states $nvc"
    nvf=`grep -i -w nvf $input | awk '{printf $2}'`
    echo "number of vf states $nvf"
    echo

    echo 'computing <0|vc>'
    cat $bendinit_pot > bpot.inp
    cat $benddecay_pot >> bpot.inp

    cat > fc_0vc.inp <<EOF
*** eSPec input file ***
========================
**MAIN
*TITLE
 +++++ cs2 Raman fin +++++
*DIMENSION
.1D
 $bend_np/
*POTENTIAL
.FILE
bpot.inp
*MASS
$bend_mass/
*TPCALC
.SPECTRUM
.TI
*INIEIGVC
.CALC
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

**TI
*TPDIAG
.MTRXDIAG 
*NIST
 1  0/
*NFST
 $nvc  0/
*ABSTOL
 1D-6/

**END

EOF

    cp fc_0vc.inp input.spc
    time $espec > fc_0vc.out

    sed -n "/Spectrum/,/The/p" fc_0vc.out | sed "/Spec/ d" | sed "/==/ d" | sed "/*/ d" | sed "/The/ d" | awk '{printf $2"\t"$4" "$5" "$6"\n"}' > fc_0vc.dat

    echo 'done!'
    echo

    echo 'computing <vc|vf>'

    cat $benddecay_pot > bpot.inp
    cat $bendfin_pot >> bpot.inp

    cat > fc_vcvf.inp <<EOF
*** eSPec input file ***
========================
**MAIN
*TITLE
 +++++ cs2 Raman fin +++++
*DIMENSION
.1D
 $bend_np/
*POTENTIAL
.FILE
bpot.inp
*MASS
$bend_mass/
*TPCALC
.SPECTRUM
.TI
*INIEIGVC
.CALC
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

**TI
*TPDIAG
.MTRXDIAG 
*NIST
 $nvc  0/
*NFST
 $nvf  0/
*ABSTOL
 1D-6/

**END

EOF

    cp fc_vcvf.inp input.spc
    time $espec > fc_vcvf.out

    sed -n "/Spectrum/,/The/p" fc_vcvf.out | sed "/Spec/ d" | sed "/==/ d" | sed "/*/ d" | sed "/The/ d" | awk '{printf $2"\t"$4" "$5" "$6"\n"}' > fc_vcvf.dat

    echo 'done!'
    echo
    echo 'Finished Franck-Condon section'
    echo

    if [ "$print_level" == "minimal" ]; then 
	rm input.spc fc_vcvf.inp fc_0vc.inp bpot.inp initial_spc.aux
    elif [ "$print_level" == "essential" ]; then
	rm input.spc fc_vcvf.inp fc_0vc.inp bpot.inp initial_spc.aux
    elif [ "$print_level" == "intermediate" ]; then
	rm input.spc fc_vcvf.inp fc_0vc.inp bpot.inp initial_spc.aux
    fi


fi

#------ Multidimensional FC factors---------------------------------------------#
if [ "$runtype" == "-all" ] || [ "$runtype" == "-fc" ] || [ "$runtype" == "-xas" ] && [ $tpmodel -eq 2 ] && [ "$dim_fc" != ".1D" ]; then
    echo
    echo "Multidimensional FC calculation"
    echo "-------------------------------"
    echo
    nmodes=$(echo $dim_fc | sed 's/[^0-9]*//g')
    for (( i=0; i<nmodes; i++ ))
    do
	# potential files
	cc=$(echo $i + 2 | bc)
	fc_mass[$i]=`grep -i -w fc_mass $input | awk -v col=$cc '{printf $col}'`
	fc_nvc[$i]=`grep -i -w fc_nvc $input | awk -v col=$cc '{printf $col}'`
	fc_nvf[$i]=`grep -i -w fc_nvf $input | awk -v col=$cc '{printf $col}'`

	fc_init_pot[$i]=`grep -i -w fc_init_pot $input | awk -v col=$cc '{printf $col}'`
	fc_decay_pot[$i]=`grep -i -w fc_decay_pot $input | awk -v col=$cc '{printf $col}'`
	fc_fin_pot[$i]=`grep -i -w fc_fin_pot $input | awk -v col=$cc '{printf $col}'`

	echo "Mode $i "
	echo "np : ${fc_np[$i]}"
	echo "mass : ${fc_mass[$i]} amu"
	echo "nvc : ${fc_nvc[$i]}"
	echo "nvf : ${fc_nvf[$i]}"
	echo "potentials:  ${fc_init_pot[$i]}, ${fc_decay_pot[$i]}, ${fc_fin_pot[$i]}"
	echo
    done


    echo 'computing Franck-Condon amplitudes'
    python -c "import functions as f; f.get_multd_fc($nmodes,'$input')"

    cat evc_temp > fc_0vc.out
    cat fc_0vc_temp >> fc_0vc.out

    sed -n "/Spectrum/,/The/p" fc_0vc.out | sed "/Spec/ d" | sed "/==/ d" | sed "/*/ d" | sed "/The/ d" | awk '{printf $2"\t"$4" "$5" "$6"\n"}' > fc_0vc.dat

    cat evf_temp > fc_vcvf.out
    cat fc_vcvf_temp >> fc_vcvf.out

    sed -n "/Spectrum/,/The/p" fc_vcvf.out | sed "/Spec/ d" | sed "/==/ d" | sed "/*/ d" | sed "/The/ d" | awk '{printf $2"\t"$4" "$5" "$6"\n"}' > fc_vcvf.dat


    echo 'done!'
    echo
    echo 'Finished Franck-Condon section'
    echo
    
    if [ "$print_level" == "minimal" ]; then 
	rm input.spc fc_vcvf.inp fc_0vc.inp bpot.inp initial_spc.aux
    elif [ "$print_level" == "essential" ]; then
	rm input.spc fc_vcvf.inp fc_0vc.inp bpot.inp initial_spc.aux
    elif [ "$print_level" == "intermediate" ]; then
	rm input.spc fc_vcvf.inp fc_0vc.inp bpot.inp initial_spc.aux
    fi   
fi


if [ "$runtype" == "-all" ] || [ "$runtype" == "-cond" ] || [ "$runtype" == "-cfin" ] ; then
    #---------------|Phi(0)> calculation--------------#
    echo 'Generating initial conditions for second propagation'
    echo

    
    if [ -f  "${jobid}.log" ]; then
	echo "# starting log file" > ${jobid}.log
    fi

    # bending states variables
    #2d+1d like run
    if [ $tpmodel -eq 0 ]; then 
	nvc=`grep -i -w nvc $input | awk '{printf $2}'`
	nvf=`grep -i -w nvf $input | awk '{printf $2}'`
    #1d or 2D run 
    elif [ $tpmodel -eq 1 ]; then 
	nvc=1
	nvf=1
    elif [ $tpmodel -eq 2 ]; then 
	nmodes=$(echo $dim_fc | sed 's/[^0-9]*//g')
	nvc=1
	nvf=1
	for (( i=0; i<nmodes; i++ ))
	do
	    # potential files
	    cc=$(echo $i + 2 | bc)
	    fc_nvc_wk=`grep -i -w fc_nvc $input | awk -v col=$cc '{printf $col}'`
	    fc_nvf_wk=`grep -i -w fc_nvf $input | awk -v col=$cc '{printf $col}'`
	    nvc=$(echo "$nvc * $fc_nvc_wk" | bc)
	    nvf=$(echo "$nvf * $fc_nvf_wk" | bc)
	done    
	echo '(FC) total number of intermediate vibrational states: ' $nvc
	echo '(FC) total number of final vibrational states: ' $nvf
    fi



    Vg_min=`grep -i Vg_min $input | awk '{printf $2}'`
    Vd_min=`grep -i Vd_min $input | awk '{printf $2}'`
    Vf_min=`grep -i Vf_min $input | awk '{printf $2}'`
    Vd=`grep -i Vd_vert $input | awk '{printf $2}'`

    if [ -f "${jobid}_init.out" ]; then

	E0=`grep '|     0        |' ${jobid}_init.out | awk '{printf $4}'`
	if [ -z "$E0" ]; then
	    E0=`grep -i -w "E0" $input | awk '{printf $2}'` 
	fi

	if [ -z "$E0" ]; then
	    echo "could not find the value of E0"
	    echo "if you are reading a initial wavefunction you did not provide E0"
	    echo "please check your input"
	    exit 666  
	else
	    echo "Initial energy" $E0
	    echo "Initial energy" $E0 >> $jobid.log
	fi
    else

	echo "failed to find initial propagation output file ${jobid}_init.out"
	echo "Attempting to read E0 from input file"
	E0=`grep -i -w "E0" $input | awk '{printf $2}'`
	if [ -z "$E0" ]; then
	    echo "could not find the value of E0"
	    echo "please check your input"
	    exit 666
	fi

    fi

    
    #------- bending |0> energy
    #2d+1d like run
    if [ $tpmodel -eq 0 ]; then 

	bE0=`sed -n "/the initial state/,/End of file/p" fc_0vc.out | grep "|     0        |" | awk '{printf $4}'`
	echo ""
	echo "Bending ground state energy bE0 = $bE0"
	check=`grep -i -w  bE0 ${jobid}.log | awk '{printf $1}'`
	if [ -z "$check" ]; then
	    echo "bE0 $bE0" >> ${jobid}.log
	fi
    #1d+nd
    elif [ $tpmodel -eq 2 ]; then 

	bE0=`sed -n "/the initial state/,/End of file/p" fc_0vc.out | grep "|     0        |" | awk '{printf $4}'`
	echo ""
	echo "mult-mode vibrational ground state energy bE0 = $bE0"
	check=`grep -i -w  bE0 ${jobid}.log | awk '{printf $1}'`
	if [ -z "$check" ]; then
	    echo "bE0 $bE0" >> ${jobid}.log
	fi

    #1d or 2D run
    elif [ $tpmodel -eq 1 ]; then  

	bE0='0.000000000'
	check=`grep -i -w  bE0 ${jobid}.log | awk '{printf $1}'`
	if [ -z "$check" ]; then
	    echo "bE0 $bE0" >> ${jobid}.log
	fi
	
    fi

    #--------

    E0tot=$(awk "BEGIN {print  $E0 + $bE0}")
    echo "E0tot" $E0tot >> $jobid.log
    echo
    echo "total initial energy E0tot = $E0tot"

    #defining detuning (for the case that omega_incoming is given as input)
    if [ -z "$all_detunings" ]; then
	echo
	echo 'defining detuning values from the inputed incoming photon energies'
	echo "detuning = omega_incoming - V_d + E0tot"
	all_detunings=`echo "$all_omega" | awk -v Vd=$Vd -v E=$E0tot '{for(i=1; i<=NF; i++) print $i-(Vd * 27.2114)+(E*27.2114)}' | awk '{printf $1" "}'`
	echo "$all_detunings eV"
    fi

    #------- 

    check=`grep -i "resonance frequency:" ${jobid}.log | awk '{printf $1}'`
    if [ -z "$check" ]; then
	omres=$(awk "BEGIN {print $Vd - $Vg_min - $E0tot}")
	echo "resonance frequency: $omres"
	echo "resonance frequency: $omres" >> $jobid.log
    else
	omres=`grep -i "resonance frequency:" ${jobid}.log | awk '{printf $3}'`
	echo "resonance frequency: $omres"
    fi

    check=`grep -i "shifted resonance frequency:" ${jobid}.log | awk '{printf $1}'`
    if [ -z "$check" ]; then

	Eres=$(awk "BEGIN {print $Vd - $Vd_min }")
	echo "shifted resonance frequency: $Eres"
	echo "shifted reson. frequency: $Eres" >> $jobid.log
    else
	Eres=`grep -i "shifted resonance frequency:" ${jobid}.log | awk '{printf $4}'`
	echo "shifted resonance frequency: $Eres"
    fi

    check=`grep -i "ground to final gap:" ${jobid}.log | awk '{printf $1}'`
    if [ -z "$check" ]; then
	Vgf_gap=$(awk "BEGIN {print $Vf_min - Vg_min}")
	echo "ground to final gap: $Vgf_gap"
	echo "ground to final gap: $Vgf_gap" >> $jobid.log
    else
	Vgf_gap=`grep -i "ground to final gap:" ${jobid}.log | awk '{printf $5}'`
	echo "ground to final gap: $Vgf_gap"	
    fi

    nfiles=`ls wf_data/ReIm_*.dat | awk '{printf $1"\n"}' | tail -1 | cut -c 14-17`
    last_file=`ls wf_data/ReIm_*.dat | awk '{printf $1"\n"}' | tail -1`
    rtime=`head -1 $last_file | awk '{printf $3}'`

    echo
    echo "number of wavepacket files: $nfiles"
    #ndetun=`echo $all_detunings | wc -w`

    #-----------------------------
    if [ -f "inwf_"*".dat" ];
    then
	echo 'removing previous initial wavepacket file!'
	rm inwf_*.dat 
    fi
    #-----------------------------
    if [ -f "intens_"*".dat" ];
    then
	echo 'removing previous wavepacket norm file (intens)!'
	rm intens_*.dat
    fi
   
 ## -------------- Loop for bending modes ----------------   
    for (( i=0 ; i < ${nvc} ; i++ )); do

	#2d+1d like run
	if [ $tpmodel -eq 0 ] || [ $tpmodel -eq 2 ]; then 
	    Evc[$i]=`sed -n "/from final state/,/Spectrum/p" fc_0vc.out | grep "|     $i        |" | awk '{printf $4}'`
	    shift[$i]=$(awk "BEGIN {print -${Evc[$i]} }")
	    echo
	    echo "bending state vc = $i, Evc$i = ${Evc[$i]} a.u."
	    echo "associated shift, -Evc$i = ${shift[$i]}"
	    echo	
	#1d or 2D run
	elif [ $tpmodel -eq 1 ]; then
	    Evc[$i]='0.0000000'
	    shift[$i]=$(awk "BEGIN {print -${Evc[$i]} }")
	    echo
	    echo "associated shift, -Evc$i = ${shift[$i]}"
	    echo	
	fi

	#transition dipole moments
	if [ "$dip_df" == "y" ]; then
	    echo ">>" $dip_df
	    echo "transition dipole moments will be read from file: $dipole_df"
	    dipole_df_input="Qoperator ../$dipole_df"
	else
	    dipole_df_input=""
	fi
	
	cat > raman.inp <<EOF
# eSPec-Raman input

*Main
dimension
$dim
$npoints
mass: $mass
filename: ../wf_data/ReIm_
nfiles: 1 $nfiles
timeinterval: 0.000 $rtime

*Propagation
Ereso: $Eres
shift: ${shift[$i]}
detuning: $all_detunings
Fourier: 10
Window
.EXPDEC $Gamma
$dipole_df_input

EOF


    # Run Raman in parallel
    rdir=${jobid}-cond-vc${i}
    if [ -d "$rdir" ]; then
        rm -r $rdir
    fi
    mkdir $rdir
    cp raman.inp $rdir 
    cd $rdir

    #2d+1d like run
    if [ $tpmodel -eq 0 ] || [ $tpmodel -eq 2 ] ; then
	echo "Running parallel |Phi(0)> calculation for Bending $i"
    #1d or 2D run
    elif [ $tpmodel -eq 1 ]; then 
	echo "Running |Phi(0)> calculation"
    fi

    time $raman > ${jobid}_raman_Evc$i.out &
    cd ../ 

    done # End parallel loop
    wait

    for (( i=0 ; i < ${nvc} ; i++ )); do
	
	rdir=${jobid}-cond-vc${i}

        # this had a problem with ${jobid} with _ in the name
	#all_detunings_files=`ls $rdir/wp_* | cat | cut  -f2 -d"_" | sed "s/.dat//" | awk '{printf $1" "}'`
	all_detunings_files=`grep -i -w detuning: $rdir/${jobid}_raman_Evc$i.out | sed 's/detuning://'`
	if [ "$i" -eq "0" ]; then
	    #------------------->>
	    check=`grep -i -w detun ${jobid}.log | awk '{printf $1}'`
	    if [ -z "$check" ]; then
		echo "detun" $all_detunings_files >> ${jobid}.log
	    fi
	    #------------------->>
	fi

	#echo "$all_detunings_files"
	# rename files according to respective bending state -----
	for detun in `echo $all_detunings_files`
	do
	    file=$rdir/wp_${detun}.dat
	    mv $file wp-vc${i}_${detun}.dat
	    cat wp-vc${i}_${detun}.dat >> inwf_${detun}.dat
	    echo " " >> inwf_${detun}.dat
	    echo "vc${i}" `head -1 wp-vc${i}_${detun}.dat | awk '{printf $8}'` >> intens_${detun}.dat
	done # End of detuning loop

    done # End of bending modes loop

    echo
    echo 'initial conditions generated!'
    echo

    # cleaning up ------
    if [ "$print_level" == "minimal" ]; then 
	rm raman.inp ${jobid}_raman_Evc*.out wp2_*
	rm -r wf_data fft_check_in.dat fft_check_out.dat
	rm ${jobid}_init.out
    elif [ "$print_level" == "essential" ]; then
	rm raman.inp ${jobid}_raman_Evc*.out wp2_*
	rm fft_check_in.dat fft_check_out.dat
    elif [ "$print_level" == "intermediate" ]; then
	rm raman.inp wp2_*
	rm fft_check_in.dat fft_check_out.dat
    elif [ "$print_level" == "savedisk" ]; then
	rm -r wf_data
	for (( i=0 ; i < ${nvc} ; i++ )); do
	    rdir=${jobid}-cond-vc${i}
	    rm $rdir/fft_spline.bcoef
	done
    fi

    #-------------------------------------------------#
fi

if [ "$runtype" == "-all" ] || [ "$runtype" == "-fin" ] || [ "$runtype" == "-cfin" ]; then

    #---------------Final Propagation-----------------#
    echo ' Starting final propagation'

    if [ -f "fcorrel.dat" ]; then
	rm fcorrel.dat
    fi

   
    # bending states variables
    #2d+1d like run
    if [ $tpmodel -eq 0 ]; then 
	nvc=`grep -i -w nvc $input | awk '{printf $2}'`
	nvf=`grep -i -w nvf $input | awk '{printf $2}'`
    #1d or 2D run
    elif [ $tpmodel -eq 1 ]; then 
	nvc=1
	nvf=1
    elif [ $tpmodel -eq 2 ]; then 
	nvc=1
	nvf=1
    fi



    work=`grep -i -w detun ${jobid}.log | awk '{printf $1}'`
    all_detunings_files=`grep -i -w detun ${jobid}.log | sed "s/\<$work\>//g"`

    for detun in `echo $all_detunings_files`
    do
	if [ -f "fcorrel_$detun.dat" ]; then
	    rm fcorrel_$detun.dat
	fi

	for ((i=0 ; i < $nvc ; i++));  do  

	    echo "running vc $i and detuning = $detun"

	    file=wp-vc${i}_${detun}.dat
	    fileinp=wp-vc${i}_${detun}.inp
	    cp $file $fileinp
	    cat $final_pot >> $fileinp

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
$fileinp
*MASS
$mass/
*TPCALC
.PROPAGATION
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
.YES
*PRTONLYREIM
$crosskey
$cross/

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
*FOURIER
 14/
$abs
$abs1
$abstren
$absrang

**END
EOF

#----------------- To run in parallel -------------------
        
    # Copying files to separate directory
        rdir=fin_vc${i}_$detun
        mkdir $rdir
        cp $fileinp $rdir
        cp input.spc $rdir
        cd $rdir

    # Running
        echo "Running vc $i in background"
	    time $espec > ${jobid}_vc${i}_$detun.out &
        cd ../
    done
    wait
    echo "All finished!"

    echo
    echo "computing correlation functions"
    echo

    for ((i=0 ; i < $nvc ; i++));  do   

        cd fin_vc${i}_$detun 

        sed "/#/ d" ../inwf_${detun}.dat > inwf.dat
	
	nfiles=`ls ReIm_*.dat | awk '{printf $1"\n"}' | tail -1 | cut -c 6-9`
	last_file=`ls ReIm_*.dat | awk '{printf $1"\n"}' | tail -1`
	rtime=`head -1 $last_file | awk '{printf $3}'`



	cat > correl.inp <<EOF
# correl input

*Main
runtype: correl
dimension
$dim
$npoints
mass: $mass
filename: ReIm_
nfiles: 1 $nfiles
timeinterval: 0.000 $rtime

*Correlation
wfunctions $nvc inwf.dat

EOF

        echo "Running fcorrel in Background for vc $i"
	time $fcorrel > ${jobid}-correl_vc${i}_$detun.out &
        cd ../
    done
    wait

    for ((i=0 ; i < $nvc ; i++));  do
        cd fin_vc${i}_$detun
	sed -n "/All desired/,/End/p" ${jobid}-correl_vc${i}_$detun.out | sed "/#/ d" | sed "/---/ d" | sed "/All/ d" > fcorrel_vc${i}_$detun.dat
        cp fcorrel_vc${i}_$detun.dat ../
	cat fcorrel_vc${i}_$detun.dat >> ../fcorrel_$detun.dat
        corr_np=`cat fcorrel_vc${i}_$detun.dat | sed '/^\s*$/d' | cat -n | tail -1 | awk '{printf $1}'`
        #corr_np=`cat -n  fcorrel_vc${i}_$detun.dat | tail -2 | awk '{printf $1" "}' | tail -2 | awk '{printf $1}'`
        cd ../
    done

# End of Detuning loop
    done

    check=`grep -i -w "corr_np" ${jobid}.log | awk '{printf $1}'`
    if [ -z "$check" ]; then
	echo "corr_np $corr_np" >> ${jobid}.log
    fi

    #rm ReIm_*

    echo "All correlation functions computed!"

    #--- cleaning up
    if [ "$print_level" == "minimal" ]; then 
	rm correl.inp  ${jobid}_vc*_*.out inwf_*.dat wp-vc*_*.inp
	rm  input.spc wp-vc*_*.dat inwf.dat ${jobid}-correl_vc*_*.out
    elif [ "$print_level" == "essential" ]; then
	rm correl.inp input.spc inwf_*.dat wp-vc*_*.inp inwf.dat
	rm ${jobid}-correl_vc*_*.out
    elif [ "$print_level" == "intermediate" ]; then
	rm correl.inp input.spc inwf_*.dat wp-vc*_*.inp inwf.dat
    elif [ "$print_level" == "savedisk" ]; then
	for detun in `echo $all_detunings_files`; do
	    for ((i=0 ; i < $nvc ; i++));  do
		rdir=fin_vc${i}_$detun
		rm $rdir/ReIm_*
	    done
	done
    fi


    #num=`cat -n  bkp-${jobid}_$detun.spec | tail -1 | awk '{printf $1}'`
    #norm=`head -1 $file | awk '{printf $8}'`
    #omres=`grep "resonance frequency:" $jobid.log | awk '{printf $3}'`
    #omres=$(awk "BEGIN {print $omres * 27.2114}")
    #Vgf_gap=`grep "ground to final gap:" $jobid.log | awk '{printf $5}'`
    #E0=`grep "Initial energy" $jobid.log | awk '{printf $3}'`
    #Vgf_gap=$(awk "BEGIN {print ($Vgf_gap - $E0) * 27.2114}")

    #-------------------------------------------------#

fi

if [ "$runtype" == "-all" ] || [ "$runtype" == "-cross" ]; then
    echo "Final cross section"

    
    # bending states variables
    #2d+1d like run
    if [ $tpmodel -eq 0 ]; then 
	nvc=`grep -i -w nvc $input | awk '{printf $2}'`
	nvf=`grep -i -w nvf $input | awk '{printf $2}'`
    #1d or 2D run
    elif [ $tpmodel -eq 1 ]; then 
	nvc=1
	nvf=1
    elif [ $tpmodel -eq 2 ]; then 
	nvc=1
	nvf=1
    fi

    work=`grep -i -w detun ${jobid}.log | awk '{printf $1}'`
    all_detunings_files=`grep -i -w detun ${jobid}.log | sed "s/\<$work\>//g"`
    corr_np=`grep -i -w corr_np ${jobid}.log | awk '{printf $2}'`


    if [ -f Evf.dat ]; then
	rm Evf.dat
    fi


    for ((i=0 ; i < $nvf ; i++)); do
	#2d+1d like run
	if [ $tpmodel -eq 0 ]; then 
	    Evf=`sed -n "/from final state/,/Spectrum/p" fc_vcvf.out | grep "|     $i        |" | awk '{printf $4}'`
	    echo $Evf >> Evf.dat
	#1d or 2D run
	elif [ $tpmodel -eq 1 ]; then 
	    Evf='0.000000'
	    echo $Evf >> Evf.dat
	fi
    done

    Evf=`cat Evf.dat | awk '{printf $1" "}'`
    rm Evf.dat

    if [ $tpmodel -eq 0 ]; then 
	echo "Final bending energies " $Evf
    fi

    for detun in `echo $all_detunings_files`
    do

	#------- variables for spectrum shift
	omres=`grep "resonance frequency:" $jobid.log | awk '{printf $3}'`

	Vgf_gap=`grep "ground to final gap:" $jobid.log | awk '{printf $5}'`
	E0=`grep "Initial energy" $jobid.log | awk '{printf $3}'`
	bE0=`grep -i -w  bE0 ${jobid}.log | awk '{printf $2}'`
	E0tot=$(awk "BEGIN {print  $E0 + $bE0}")

	detunau=$(awk "BEGIN {print $detun / 27.2114}")

	omega=$(awk "BEGIN {print $omres + $detunau}")
	shift=$(awk "BEGIN {print $Vgf_gap + $E0tot }")
	#------------------------------------
	#2d+1d like run
	if [ $tpmodel -eq 0 ]; then 
	    cat  fc_0vc.dat | awk '{printf $1" \n"}' > fcond.dat
	    cat  fc_vcvf.dat | awk '{printf $1" \n"}' >> fcond.dat
	    cat intens_$detun.dat | awk '{printf $2" \n"}' >> fcond.dat
	#1d or 2D run
	elif [ $tpmodel -eq 1 ]; then 
	    echo '1.00000' > fcond.dat
	    echo '1.00000' >> fcond.dat
	    cat intens_$detun.dat | awk '{printf $2" \n"}' >> fcond.dat
	fi
	   

	cat > correl.inp <<EOF
# final cross section input

*Main
runtype: spectrum
corr_np $corr_np

*crosssection
vc $nvc
vf $nvf
Evf $Evf
franckcondon fcond.dat
fcorrel fcorrel_$detun.dat
Fourier 11
omega $omega
shift $shift
Window
.SGAUSS $window

EOF

	time $fcorrel > ${jobid}-final_csection_$detun.out

	#sed -n "/Final spectrum/,/End/p"  ${jobid}-final_csection_$detun.out | sed "/#/ d" | sed "/--/ d" | sed "/Final/ d" | sed "/(eV)/ d" | sed '/^\s*$/d' > temp

	echo "# spectrum as function of emitted photon energy E', omega= $omega " > ${jobid}_$detun.spec
	sed -n "/> RIXS spectrum as function of the emitted photon energy /,/> RIXS spectrum as function of energy loss/p" ${jobid}-final_csection_$detun.out | sed "/#/ d" | sed "/--/ d" | sed "/RIXS/ d" | sed "/(eV)/ d" | sed '/^\s*$/d' | awk '{printf $1" "$2"\n"}' >> ${jobid}_$detun.spec


	echo "# spectrum as function of energy loss E - E', omega= $omega " > ${jobid}_${detun}_Eloss.spec
	sed -n "/> RIXS spectrum as function of energy loss/,/# End of Calculation/p" ${jobid}-final_csection_$detun.out | sed "/#/ d" | sed "/--/ d" | sed "/RIXS/ d" | sed "/(eV)/ d" | sed '/^\s*$/d' | awk '{printf $1" "$2"\n"}' >> ${jobid}_${detun}_Eloss.spec

	# if [ "$doshift"=="y" ]; then
	#     echo "shifting spectrum, omega=$omega eV"
	#     shift=`awk "BEGIN {print $omega -$Vgf_gap + $bE0}"`
	#     while read x y discard; do
	# 	w=$(awk "BEGIN {print $shift - $x}")
	# 	Int=$(awk "BEGIN {print $y}")
	# 	#Int=$(awk "BEGIN {print $norm * $y}")
	# 	echo "$w $Int" >> ${jobid}_$detun.spec
	#     done < temp
	#     #cat temp | awk '{printf $1" "$2"\n"}' > bkp-${jobid}_$detun.spec
	#     #cat temp | awk "BEGIN {printf $shift - $1}" >> ${jobid}_$detun.spec
	# else
	#     cat temp | awk '{printf $1" "$2"\n"}' > ${jobid}_$detun.spec   
	# fi
	
	#rm temp

	echo
	echo "Final spectrum saved to ${jobid}_$detun.spec and ${jobid}_${detun}_Eloss.spec"
	echo

    done

    #--- cleaning up
    if [ "$print_level" == "minimal" ]; then 
	rm correl.inp ${jobid}-final_csection_*.out fcond.dat fcorrel_*dat
	rm intens_*.dat fc_vcvf.out fc_0vc.out
    elif [ "$print_level" == "essential" ]; then
	rm correl.inp ${jobid}-final_csection_*.out fcond.dat
    elif [ "$print_level" == "intermediate" ]; then
	rm correl.inp fcond.dat
    fi

fi


### -------- 2D+1D XAS

if [ "$runtype" == "-xas" ] || [ "$runtype" == "-xascs" ]; then
    echo "XAS cross section"
    echo

    Vg_min=`grep -i Vg_min $input | awk '{printf $2}'`
    Vd_min=`grep -i Vd_min $input | awk '{printf $2}'`
    Vf_min=`grep -i Vf_min $input | awk '{printf $2}'`
    Vd=`grep -i Vd_vert $input | awk '{printf $2}'`

    #2d+1d like run
    if [ $tpmodel -eq 0 ]; then
	nvc=`grep -i -w nvc $input | awk '{printf $2}'`
	echo "number of core-excited bending modes included, nvc: $nvc"
	nvf=`grep -i -w nvf $input | awk '{printf $2}'`
    #1d or 2D run
    elif [ $tpmodel -eq 1 ]; then 	
	nvc='1'
	nvf='1'
    fi
    
    corr_np=`cat -n xas-fcorrel.dat | tail -1 | awk '{printf $1}'`

    if [ -f Evc.dat ]; then
	rm Evc.dat
    fi

    for ((i=0 ; i < $nvf ; i++)); do
	
	#2d+1d like run
	if [ $tpmodel -eq 0 ]; then
	    Evc=`sed -n "/from final state/,/Spectrum/p" fc_0vc.out | grep "|     $i        |" | awk '{printf $4}'`
	    echo $Evc >> Evc.dat
	#1d or 2D run
	elif [ $tpmodel -eq 1 ]; then 
	    Evc='0.000000'
	    echo $Evc >> Evc.dat
	fi
	
    done

    Evc=`cat Evc.dat | awk '{printf $1" "}'`
    rm Evc.dat

    if [ $tpmodel -eq 0 ]; then
	echo "Core-excited bending energies " $Evc "(a.u.)"
    fi

#---------
    check=`grep -i "shifted resonance frequency:" ${jobid}.log | awk '{printf $1}'`
    if [ -z "$check" ]; then
	Eres=$(awk "BEGIN {print $Vd - $Vd_min }")
	echo "shifted resonance frequency(Delta): $Eres a.u."
	echo "shifted reson. frequency: $Eres" >> $jobid.log
    else
	Eres=`grep -i "shifted resonance frequency:" ${jobid}.log | awk '{printf $4}'`
	echo "shifted resonance frequency (Delta): $Eres"
    fi
    #-------------------
    Vg_min=`grep -i Vg_min $input | awk '{printf $2}'`
    Vd=`grep -i Vd_vert $input | awk '{printf $2}'`

    #----------------------
    #    Energy of initial state
    #---------------------
    if [ -f "${jobid}_init.out" ]; then

	E0=`grep '|     0        |' ${jobid}_init.out | awk '{printf $4}'`
	if [ -z "$E0" ]; then
	    E0=`grep -i -w "E0" $input | awk '{printf $2}'` 
	fi

	if [ -z "$E0" ]; then
	    echo "could not find the value of E0"
	    echo "if you are reading a initial wavefunction you did not provide E0"
	    echo "please check your input"
	    exit 666  
	else
	    echo "Initial stretching energy $E0 a.u."
	fi
    else

	echo "failed to find initial propagation output file ${jobid}_init.out"
	echo "Attempting to read E0 from input file"
	E0=`grep -i -w "E0" $input | awk '{printf $2}'`
	if [ -z "$E0" ]; then
	    echo "could not find the value of E0"
	    echo "please check your input"
	    exit 666
	fi

    fi
    #------- bending |0> energy
    #2d+1d like run
    if [ $tpmodel -eq 0 ]; then 

	bE0=`sed -n "/the initial state/,/End of file/p" fc_0vc.out | grep "|     0        |" | awk '{printf $4}'`
	echo
	echo "Bending ground state energy bE0 = $bE0 a.u."
	check=`grep -i -w  bE0 ${jobid}.log | awk '{printf $1}'`
	if [ -z "$check" ]; then
	    echo "bE0 $bE0" >> ${jobid}.log
	fi

    #1d or 2D run
    elif [ $tpmodel -eq 1 ]; then  
	
	bE0='0.000000000'
	check=`grep -i -w  bE0 ${jobid}.log | awk '{printf $1}'`
	if [ -z "$check" ]; then
	    echo "bE0 $bE0" >> ${jobid}.log
	fi

    fi


    #--------
    E0tot=$(awk "BEGIN {print  $E0 + $bE0}")
    echo "Total initial energy E0tot = $E0tot a.u."
    #-------
    shift=$(awk "BEGIN {print $Vd - $Vg_min -$E0tot}")
    echo "XAS resonance frequency: $shift a.u."

#-------------------
    if [ -f xas-fcorrel.dat ]; then
	echo "XAS correlation function found: xas-fcorrel.dat"
    else
	echo "XAS correlation function file (xas-fcorrel.dat) not found!!"
	echo "have you run the initial propagation step?"
    fi
#-------------------
    echo
    echo
    echo "starting final spectrum calculation"
    echo

    
    #2d+1d like run
    if [ $tpmodel -eq 0 ]; then 
	cat  fc_0vc.dat | awk '{printf $1" \n"}' > fcond.dat #FOR XAS WE ONLY NEED GROUND->CORE FC FACTORS
    #1d or 2D run
    elif [ $tpmodel -eq 1 ]; then 
	 echo '1.00000' > fcond.dat
    fi
	
    window=$(awk "BEGIN {print $Gamma / 27.2114}")

	cat > correl.inp <<EOF
# XAS cross section input

*Main
runtype: xas-spectrum
corr_np $corr_np

*crosssection
vc $nvc
Evc $Evc
franckcondon fcond.dat
fcorrel xas-fcorrel.dat
Delta $Eres
shift $shift
Fourier 11
Window
.EXPDEC $window

EOF

	#time $fcorrel
	#echo
	#echo
	#echo '-----------------'

	time $fcorrel > ${jobid}-xas_csection_$detun.out

	echo "# XAS spectrum as function of photon energy" > ${jobid}_xas.spec
	sed -n "/> XAS spectrum as function od photon/,/End/p"  ${jobid}-xas_csection_$detun.out | sed "/#/ d" | sed "/--/ d" | sed "/Final/ d" | sed "/(eV)/ d" | sed '/^\s*$/d' | sed "/XAS/ d" | awk '{printf $1" "$2"\n"}' >> ${jobid}_xas.spec

	echo "# XAS spectrum as function of detuning" > ${jobid}_xas-detuning.spec
	sed -n "/> XAS spectrum as function of detuning/,/> XAS spectrum as function od photon/p"  ${jobid}-xas_csection_$detun.out | sed "/#/ d" | sed "/--/ d" | sed "/Final/ d" | sed "/(eV)/ d" | sed '/^\s*$/d' | sed "/XAS/ d" | awk '{printf $1" "$2"\n"}' >> ${jobid}_xas-detuning.spec


	
	#rm temp

	echo
	echo "XAS spectrum saved to ${jobid}_xas.spec and ${jobid}_xas-detuning.spec"
	echo


    #--- cleaning up
    if [ "$print_level" == "minimal" ]; then 
	rm correl.inp ${jobid}-xas_csection_*.out fcond.dat fcorrel_*dat
	rm intens_*.dat fc_vcvf.out fc_0vc.out
    elif [ "$print_level" == "essential" ]; then
	rm correl.inp ${jobid}-xas_csection_*.out fcond.dat
    elif [ "$print_level" == "intermediate" ]; then
	rm correl.inp fcond.dat
    fi

fi

#
#----------------- self-absorption part
#

if [ "$runtype" == "-self" ]; then
    echo "REXS cross section with Self-Absorption factor"
    echo

    if [ -f ${jobid}_xas.spec ]; then
	echo "XAS spectrum file found: ${jobid}_xas.spec"
	echo
	cat ${jobid}_xas.spec | sed "/#/ d" > temp_xas.spec
	nxas=`cat -n temp_xas.spec | tail -1 | awk '{printf $1}'`
    else
	echo "ERROR!!"	
	echo "XAS spectrum file not found!!"
	echo "Have you run a -xas calculation previously?"
	exit
    fi
#
#-----------#
#
# I need a loop on all detuning files!!!!
#    try to use the Eloss files
#    all_detunings_files=`ls ${jobid}_*.spec | cat | cut  -f2 -d"_" | sed "s/.spec//" | awk '{printf $1" "}'`
    work=`grep -i -w detun ${jobid}.log | awk '{printf $1}'`
    all_detunings_files=`grep -i -w detun ${jobid}.log | sed "s/\<$work\>//g"`

    for detun in  `echo $all_detunings_files`
    do
	#------- variables for spectrum shift
	omres=`grep "resonance frequency:" $jobid.log | awk '{printf $3}'`

	Vgf_gap=`grep "ground to final gap:" $jobid.log | awk '{printf $5}'`
	E0=`grep "Initial energy" $jobid.log | awk '{printf $3}'`
	bE0=`grep -i -w  bE0 ${jobid}.log | awk '{printf $2}'`
	E0tot=$(awk "BEGIN {print  $E0 + $bE0}")

	detunau=$(awk "BEGIN {print $detun / 27.2114}")

	omega=$(awk "BEGIN {print $omres + $detunau}")
	shift=$(awk "BEGIN {print $Vgf_gap + $E0tot }")
	#------------------------------------

	if [ -f ${jobid}_$detun.spec ]; then
	    echo
	    echo "REXS spectrum file for detuning = $detun found: ${jobid}_$detun.spec"
	    echo
	    cat ${jobid}_$detun.spec | sed "/#/ d" > temp_rexs.spec
	    nrexs=`cat -n temp_rexs.spec | tail -1 | awk '{printf $1}'`
	else
	    echo
	    echo "ERROR!!"	
	    echo "REXS spectrum file for detuning = $detun not found!!"
	    echo "Have you run a -all calculation previously?"
	    exit
	fi

	
	cat > correl.inp <<EOF
# XAS cross section input

*Main
runtype: self-abs

*crosssection
rexs-cs $nrexs temp_rexs.spec 
xas-cs  $nxas  temp_xas.spec
omega $omega

EOF

	time $fcorrel > ${jobid}-rexs-sa_csection_$detun.out
	echo "# spectrum as function of emitted photon energy E', omega= $omega " > ${jobid}_$detun-sa.spec
	sed -n "/> REXS-SA spectrum as function of emitted photon energy/,/> REXS-SA spectrum as function of energy loss/p" ${jobid}-rexs-sa_csection_$detun.out | sed "/#/ d" | sed "/--/ d" | sed "/REXS/ d" | sed "/(eV)/ d" | sed '/^\s*$/d' | awk '{printf $1" "$2"\n"}' >> ${jobid}_$detun-sa.spec


	echo "# spectrum as function of energy loss E - E', omega= $omega " > ${jobid}_${detun}_Eloss-sa.spec
	sed -n "/> REXS-SA spectrum as function of energy loss/,/# End of Calculation/p" ${jobid}-rexs-sa_csection_$detun.out | sed "/#/ d" | sed "/--/ d" | sed "/REXS/ d" | sed "/(eV)/ d" | sed '/^\s*$/d' | awk '{printf $1" "$2"\n"}' >> ${jobid}_${detun}_Eloss-sa.spec
	#---------
    done
    #--------------------
fi


echo
echo 'eSPec-Raman script finished'
echo
date
echo
