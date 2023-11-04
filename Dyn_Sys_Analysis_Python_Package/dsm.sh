#! /bin/bash
# dir: data directory.
# xf: x data file name.
# yf: y data file name.
# ouname: ouput file name.
# xvar: x variable name.
# yvar: y variable name.

i=$1
export dir=""
export pythondir=""

if [ i -eq 2 ]; then

	export xf=""
	export yf=""
	export xvar=""
	export yvar=""
	export ouname=""

	${pythondir}  Main_Bivariate_Analysis.py
	unset xf;unset yf;unset xvar;unset yvar;unset ouname
else
	export xf=""
	export xvar=""
	export ouname=""

	${pythondir}  Main_Univariate_Analysis.py

	unset xf;unset xvar
fi

unset dir;unset pythondir



