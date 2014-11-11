#!/bin/csh -f

# Defining colours for shell
set GREEN  = "\033[1;32m"
set RED    = "\033[1;31m"
set PINK   = "\033[1;35m"
set BLUE   = "\033[1;34m"
set YELLOW = "\033[1;33m"
set CYAN   = "\033[1;36m"
set NORMAL = "\033[0;39m"

# Configuring MA5 environment variable
setenv MA5_BASE /Users/lb8075/madanalysis5

# Configuring PATH environment variable
if ( $?PATH ) then
setenv PATH $MA5_BASE/tools/fastjet/bin:"$PATH":/Users/lb8075/SecondRoot/root/bin
else
setenv PATH $MA5_BASE/tools/fastjet/bin:/Users/lb8075/SecondRoot/root/bin
endif

# Configuring LD_LIBRARY_PATH environment variable
if ( $?LD_LIBRARY_PATH ) then
setenv LD_LIBRARY_PATH $MA5_BASE/tools/SampleAnalyzer/Lib/:$MA5_BASE/tools/fastjet/lib:"$LD_LIBRARY_PATH":/Users/lb8075/SecondRoot/root/lib:/usr/lib
else
setenv LD_LIBRARY_PATH $MA5_BASE/tools/SampleAnalyzer/Lib/:$MA5_BASE/tools/fastjet/lib:/Users/lb8075/SecondRoot/root/lib:/usr/lib
endif

# Configuring DYLD_LIBRARY_PATH environment variable
if ( $?DYLD_LIBRARY_PATH ) then
setenv DYLD_LIBRARY_PATH $MA5_BASE/tools/SampleAnalyzer/Lib/:$MA5_BASE/tools/fastjet/lib:"$DYLD_LIBRARY_PATH":/Users/lb8075/SecondRoot/root/lib:/usr/lib
else
setenv DYLD_LIBRARY_PATH $MA5_BASE/tools/SampleAnalyzer/Lib/:$MA5_BASE/tools/fastjet/lib:/Users/lb8075/SecondRoot/root/lib:/usr/lib
endif

# Checking that all environment variables are defined
if ( $?MA5_BASE && $?PATH && $?LD_LIBRARY_PATH && $?DYLD_LIBRARY_PATH ) then
echo $YELLOW"--------------------------------------------------------"
echo "    Your environment is properly configured for MA5     "
echo "--------------------------------------------------------"$NORMAL
endif
