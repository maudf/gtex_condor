#!/bin/bash

### Set variables

### START OF CODE GENERATED BY Argbash v2.9.0 one line above ###
# Argbash is a bash code generator used to get arguments parsing right.
# Argbash is FREE SOFTWARE, see https://argbash.io for more info
# Generated online by https://argbash.io/generate


die()
{
	local _ret="${2:-1}"
	test "${_PRINT_HELP:-no}" = yes && print_help >&2
	echo "$1" >&2
	exit "${_ret}"
}


begins_with_short__arg_option()
{
	local first__arg_option all_short__arg_options='posfdh'
	first__arg_option="${1:0:1}"
	test "$all_short__arg_options" = "${all_short__arg_options/$first__arg_option/}" && return 1 || return 0
}

# THE DEFAULTS INITIALIZATION - OPTIONALS
_arg_path=~/Documents/Softwares/SelectionHapStats/
_arg_optimum=10
_arg_selectcoeff=0.78
_arg_frequency=0.05
_arg_directory=./


print_help()
{
	printf '%s\n' "The general script's help msg"
	printf 'Usage: %s [-p|--path <arg>] [-o|--optimum <arg>] [-s|--selectcoeff <arg>] [-f|--frequency <arg>] [-d|--directory <arg>] [-h|--help]\n' "$0"
	printf '\t%s\n' "-o, --path: Value of _arg_path path to SelectionHapStats software (default: ~/Documents/Softwares/SelectionHapStats/)"
	printf '\t%s\n' "-o, --optimum: Value of _arg_optimum phenotype (default: 10)"
	printf '\t%s\n' "-s, --selectcoeff: Selection coefficient of QTLs (default: 0.78)"
	printf '\t%s\n' "-f, --frequency: Frequency of QTLs at selection start (default: 0.005)"
	printf '\t%s\n' "-d, --directory: Path to results directory (default: ./)"
	printf '\t%s\n' "-h, --help: Prints help"
}


parse_commandline()
{
	while test $# -gt 0
	do
		_key="$1"
		case "$_key" in
		       	-p|--path)
				test $# -lt 2 && die "Missing value for the _arg_path argument '$_key'." 1
				_arg_path="$2"
				shift
				;;
			--path=*)
				_arg_path="${_key##--_arg_path=}"
				;;
			-p*)
				_arg_path="${_key##-p}"
				;;
			-o|--optimum)
				test $# -lt 2 && die "Missing value for the _arg_optimum argument '$_key'." 1
				_arg_optimum="$2"
				shift
				;;
			--optimum=*)
				_arg_optimum="${_key##--_arg_optimum=}"
				;;
			-o*)
				_arg_optimum="${_key##-o}"
				;;
			-s|--selectcoeff)
				test $# -lt 2 && die "Missing value for the _arg_selectcoeff argument '$_key'." 1
				_arg_selectcoeff="$2"
				shift
				;;
			--selectcoeff=*)
				_arg_selectcoeff="${_key##--selectcoeff=}"
				;;
			-s*)
				_arg_selectcoeff="${_key##-s}"
				;;
			-f|--frequency)
				test $# -lt 2 && die "Missing value for the _arg_frequency argument '$_key'." 1
				_arg_frequency="$2"
				shift
				;;
			--frequency=*)
				_arg_frequency="${_key##--frequency=}"
				;;
			-f*)
				_arg_frequency="${_key##-f}"
				;;
			-d|--directory)
				test $# -lt 2 && die "Missing value for the _arg_directory argument '$_key'." 1
				_arg_directory="$2"
				shift
				;;
			--directory=*)
				_arg_directory="${_key##--directory=}"
				;;
			-d*)
				_arg_directory="${_key##-d}"
				;;
			-h|--help)
				print_help
				exit 0
				;;
			-h*)
				print_help
				exit 0
				;;
			*)
				_PRINT_HELP=yes die "FATAL ERROR: Got an unexpected argument '$1'" 1
				;;
		esac
		shift
	done
}

parse_commandline "$@"

# OTHER STUFF GENERATED BY Argbash

### END OF CODE GENERATED BY Argbash (sortof) ### ])
# [ <-- needed because of Argbash


echo "Value of --_arg_path: $_arg_path"
echo "Value of --_arg_optimum: $_arg_optimum"
echo "Value of --selectcoeff: $_arg_selectcoeff"
echo "Value of --frequency: $_arg_frequency"
echo "Value of --directory: $_arg_directory"

# ] <-- needed because of Argbash

###Build output folders



if [[ "${_arg_directory}" =~ Neutral ]]
then
    folder="$_arg_directory"/ES_"$_arg_selectcoeff"_f_"$_arg_frequency"/
elif [[ "${_arg_directory}" =~ Directional ]]
     folder="$_arg_directory"/Opt_"$_arg_optimum"_ES_"$_arg_selectcoeff"_f_"$_arg_frequency"/
else
    folder="$_arg_directory"/Opt_"$_arg_optimum"_ES_"$_arg_selectcoeff"_f_"$_arg_frequency"/
fi

for i in $(seq 1 1 100)
do
    echo "Running Rep ${i}..."
    for j in $(seq 50001 100000 1900000)
    do
	echo -n ${i}$'\t'${j}$'\t' >> $folder/Opt_${opt}_ES_${ES}_f_${freq}/Rep_${i}_OPT_${opt}_ES_${ES}_f_${freq}_${tag}.results.h12.txt 
	python H12_H2H1.py $folder/Opt_${opt}_ES_${ES}_f_${freq}/Rep_${i}_OPT_${opt}_ES_${ES}_f_${freq}_${tag}.h12.txt 400 --singleWindow ${j} --window 200 >> $folder/Opt_${opt}_ES_${ES}_f_${freq}/Rep_${i}_OPT_${opt}_ES_${ES}_f_${freq}_${tag}.results.h12.txt 
    done
done